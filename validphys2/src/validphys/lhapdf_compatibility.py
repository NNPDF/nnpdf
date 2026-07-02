"""
Module for LHAPDF compatibility backends

If LHAPDF is installed, the module will transparently hand over everything to LHAPDF.
If LHAPDF is not available, it will try to use a combination of the packages
    `lhapdf-management` and `pdfflow`
which cover all the features of LHAPDF used during the fit (and likely most of validphys).

The NeoPDF interpolation library can be selected by setting ``pdf_backend: neopdf``
in the NNPDF profile (``nnprofile.yaml``), or via the ``NNPDF_PDF_BACKEND`` environment
variable which takes precedence over the profile.
"""

from functools import cached_property
import os

import numpy as np

try:
    import lhapdf

    USING_LHAPDF = True
except ModuleNotFoundError:
    import logging

    import lhapdf_management as lhapdf

    log = logging.getLogger(__name__)
    log.warning("LHAPDF was not found, using an alternative backend")

    USING_LHAPDF = False

_BACKEND_ENV_VAR = "NNPDF_PDF_BACKEND"
_VALID_BACKENDS = ("lhapdf", "pdfflow", "neopdf")


class InvalidPDFBackend(Exception):
    pass


def _active_backend():
    """Return the active PDF backend.

    Resolution order (highest priority first):
    1. ``NNPDF_PDF_BACKEND`` environment variable
    2. ``pdf_backend`` key in the NNPDF profile (``nnprofile.yaml``)
    """
    backend = os.environ.get(_BACKEND_ENV_VAR)
    if backend is None:
        try:
            from nnpdf_data.utils import get_nnpdf_profile

            backend = get_nnpdf_profile().get("pdf_backend", "lhapdf")
        except Exception:
            backend = "lhapdf"
    return backend.lower()


class _PDFFlowPDF:
    """Wrapper around the PDFFlow PDF so that it can be used as an LHAPDF
    set by validphys
    Takes as input a pdf_meta object (which is a PDFset from lhapdf_management
    and which knows where the PDF needs to be loaded from) and a single member

    Loading the PDF is done in a lazy manner since most of the time only a few
    members are needed.

    Since PDFFlow is only utilized to load the PDF for interpolation, the import
    is delayed until the first call to `mkPDF`. This allows the usage of most of
    validphys without tensorflow.
    """

    def __init__(self, pdf_meta, member):
        if USING_LHAPDF:
            raise ValueError("PDFFlow should not be instantiated when using LHAPDF")

        self._pdf_meta = pdf_meta
        self._m = member
        self._pdf = None
        self._flavors = self._pdf_meta.info["Flavors"]

    @cached_property
    def pdf(self):
        # Don't import PDF Flow until you really needed it
        import pdfflow

        if self._pdf is None:
            pdf_def = f"{self._pdf_meta.name}/{self._m}"
            self._pdf = pdfflow.mkPDF(pdf_def, self._pdf_meta.path.parent)
        return self._pdf

    def flavors(self):
        return self._flavors

    def _xfxQ_all_pid(self, x, q):
        x = np.atleast_1d(x)
        q = np.atleast_1d(q)

        res = self.pdf.py_xfxQ2_allpid(x, q**2).numpy()
        return dict(zip(self._flavors, res.T))

    def xfxQ(self, a, b, c=None):
        """Wrapper for the LHAPDF xfxQ function
        This is an overloaded function in LHAPDF so depending
        on the number of arguments we will do:
            xfxQ(flavours, x, Q)
        or
            xfxQ(x, q)

        All of x/q/flavours can be either a scalar or an array
        """
        if c is None:
            return self._xfxQ_all_pid(a, b)

        # PDFFlow doesn't allow to ask for flavours that do not exist
        # so let us retrieve all and return 0s for non existing flavs
        ret_dict = self.xfxQ(b, c)
        zeros = np.zeros_like(b)

        if isinstance(a, int):
            return ret_dict.get(a, zeros)
        return np.array([ret_dict.get(i, zeros) for i in a]).T

    def xfxQ2(self, a, b, c=None):
        """Wrapper for LHAPDF xfxQ2 function, like xfxQ for Q2"""
        if c is None:
            return self.xfxQ(a, np.sqrt(b))
        return self.xfxQ(a, b, np.sqrt(c))


class _NeoPDFPDF:
    """Thin wrapper around a single NeoPDF member exposing the LHAPDF-compatible interface."""

    def __init__(self, neo_member):
        self._member = neo_member
        self._pids = neo_member.pids()

    def flavors(self):
        return self._pids

    def _xfxQ_all_pid(self, x, q):
        scalar_input = np.ndim(x) == 0 and np.ndim(q) == 0
        x = np.atleast_1d(x)
        q = np.atleast_1d(q)
        vals = np.array(
            [
                self._member.xfxQ2_allpids(self._pids, float(xi), float(qi) ** 2)
                for xi, qi in zip(x, q)
            ]
        )  # (n_points, n_pids)
        if scalar_input:
            return dict(zip(self._pids, vals[0]))
        return dict(zip(self._pids, vals.T))

    def xfxQ(self, a, b, c=None):
        if c is None:
            return self._xfxQ_all_pid(a, b)
        ret_dict = self.xfxQ(b, c)
        zeros = np.zeros_like(b)
        if isinstance(a, int):
            return ret_dict.get(a, zeros)
        return np.array([ret_dict.get(i, zeros) for i in a]).T

    def xfxQ2(self, a, b, c=None):
        if c is None:
            return self.xfxQ(a, np.sqrt(b))
        return self.xfxQ(a, b, np.sqrt(c))


def make_pdf(pdf_name, member=None):
    """Load a single member if specified, otherwise load the entire set as a list.

    If LHAPDF is provided, it returns LHAPDF PDF instances otherwise it returns and
    object which is _compatible_ with LHAPDF for lhapdf functions for the selected
    backend.

    The backend can be overridden by setting the ``NNPDF_PDF_BACKEND`` environment
    variable to ``neopdf`` to use the NeoPDF interpolation library.

    Parameters:
    -----------
        pdf_name: str
            name of the PDF to load
        member: int
            index of the member of the PDF to load

    Returns:
    --------
        list(pdf_sets)
    """
    backend = _active_backend()

    if backend not in _VALID_BACKENDS:
        raise InvalidPDFBackend(f"Unknown backend {backend!r}. Options are: {_VALID_BACKENDS}")

    if backend == "lhapdf":
        if member is None:
            return lhapdf.mkPDFs(pdf_name)
        return [lhapdf.mkPDF(pdf_name, member)]

    if backend == "neopdf":
        from neopdf.pdf import PDF as _NeoPDF

        members = _NeoPDF.mkPDFs(pdf_name)
        if member is None:
            return [_NeoPDFPDF(m) for m in members]
        return [_NeoPDFPDF(members[member])]

    # backend == "pdfflow"
    pdf_meta = lhapdf.load_pdf_meta(pdf_name)
    if member is None:
        return [_PDFFlowPDF(pdf_meta, m) for m in range(len(pdf_meta))]
    return [_PDFFlowPDF(pdf_meta, member)]
