"""
    Module for LHAPDF compatibility backends

    If LHAPDF is installed, the module will transparently hand over everything to LHAPDF
    if LHAPDF is not available, it will try to use a combination of the packages
        `lhapdf-management` and `pdfflow`
    which cover all the features of LHAPDF used during the fit (and likely most of validphys)

    Eventually this module will allow us to transition to an under-development python/rust
    PDF interpolation library.
"""
import numpy as np

try:
    import lhapdf

    USING_LHAPDF = True
except ModuleNotFoundError:
    import logging

    import lhapdf_management as lhapdf
    import pdfflow

    log = logging.getLogger(__name__)
    log.warning("LHAPDF was not found, using an alternative backend")

    USING_LHAPDF = False


class _PDFFlowPDF:
    """Wrapper around the PDFFlow PDF so that it can be used as an LHAPDF
    set by validphys
    Takes as input a pdf_meta object (which is a PDFset from lhapdf_management
    and which knows _where_ the PDF needs to be loaded from) and a single member

    Loading the PDF is done in a lazy manner since most of the time only a few members are needed.
    """

    def __init__(self, pdf_meta, member):
        if USING_LHAPDF:
            raise ValueError("PDFFlow should not be instantiated when using LHAPDF")

        self._pdf_meta = pdf_meta
        self._m = member
        self._pdf = None

    @property
    def pdf(self):
        if self._pdf is None:
            pdf_def = f"{self._pdf_meta.name}/{self._m}"
            self._pdf = pdfflow.mkPDF(pdf_def, self._pdf_meta.path.parent)
        return self._pdf

    @property
    def flavors(self):
        return self._pdf_meta.info["Flavors"]

    def _xfxQ_all_pid(self, x, q):
        if isinstance(x, float):
            x = np.array([x])
        if isinstance(q, float):
            q = np.array([q])

        res = self.pdf.py_xfxQ2_allpid(x, q**2).numpy()
        return dict(zip(self.flavors, res.T))

    def xfxQ(self, a, b, c=None):
        """Wrapper for the LHAPDF xfxQ function
        This is an overloaded function in LHAPDF so depending
        on the number of arguments we will do:
            xfxQ(flavours, x, Q)
        or
            xfxQ(x, q)

        And x/q/flavours can be either an scalar or an array
        """
        if c is None:
            return self._xfxQ_all_pid(a, b)

        # PDFFlow doesn't allow to ask for flavours that do not exist
        ret_dict = self.xfxQ(b, c)
        zeros = np.zeros_like(b)
        return [ret_dict.get(i, zeros) for i in a]


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
        return [_PDFFlowPDF(pdf_meta, m) for m in len(pdf_meta)]
    return [_PDFFlowPDF(pdf_meta, member)]
