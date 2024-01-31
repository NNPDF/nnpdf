"""
    Module for LHAPDF compatibility backends

    If LHAPDF is installed, the module will transparently hand over everything to LHAPDF
    if LHAPDF is not available, it will try to use a combination of the packages
        `lhapdf-management` and `pdfflow`
    which cover all the features of LHAPDF used during the fit (and likely most of validphys)
"""
from functools import cached_property

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


class _PDFFlowPDF:
    """Wrapper around the PDFFlow PDF so that it can be used as an LHAPDF
    set by validphys
    Takes as input a pdf_meta object (which is a PDFset from lhapdf_management
    and which knows _where_ the PDF needs to be loaded from) and a single member

    Loading the PDF is done in a lazy manner since most of the time only a few members are needed.

    Since PDFFlow is only utilized to load the PDF for interpolation, the import is delayed until
    the first call to `mkPDF`. This allows the usage of most of validphys without tensorflow.
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
        return [ret_dict.get(i, zeros) for i in a]

    def xfxQ2(self, a, b, c=None):
        """Wrapper for LHAPDF xfxQ2 function, like xfxQ for Q2"""
        if c is None:
            return self.xfxQ(a, np.sqrt(b))
        return self.xfxQ(a, b, np.sqrt(c))


def make_pdf(pdf_name, member=None):
    """Load a PDF
    if member is given, load the single member otherwise, load the entire set as a list

    if LHAPDF is provided, it returns LHAPDF PDF instances
    otherwise it returns and object which is _compatible_ with LHAPDF
    for lhapdf functions for the selected backend

    Parameters:
    ----------
        pdf_name: str
            name of the PDF to load
        member: int
            index of the member of the PDF to load

    Returns:
    -------
        list(pdf_sets)
    """
    if USING_LHAPDF:
        if member is None:
            return lhapdf.mkPDFs(pdf_name)
        return [lhapdf.mkPDF(pdf_name, member)]

    pdf_meta = lhapdf.load_pdf_meta(pdf_name)
    if member is None:
        return [_PDFFlowPDF(pdf_meta, m) for m in range(len(pdf_meta))]
    return [_PDFFlowPDF(pdf_meta, member)]
