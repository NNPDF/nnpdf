"""
    Module containing an LHAPDF class compatible with validphys
    using the official lhapdf python interface.

    It exposes an interface mostly compatible with libNNPDF::LHAPDFSet
    so it can be used as a drop-in replacement.

    The ``.members`` and ``.central_member`` of the ``LHAPDFSet`` are
    LHAPDF objects (the typical output from ``mkPDFs``) and can be used normally.
    For MC PDFs the ``central_member`` is the average of all replicas (members 1-100)
    while for Hessian PDFfs the ``central_member`` is also ``members[0]``

    Examples
    --------
    >>> from validphys.lhapdfset import LHAPDFSet
    >>> pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", "replicas")
    >>> len(pdf.members)
    100
    >>> pdf.central_member.alphasQ(91.19)
    0.11800
    >>> pdf.members[0].xfxQ2(0.5, 15625)
    {-5: 6.983360500601136e-05,
    -4: 0.0021818063617227604,
    -3: 0.00172453472243952,
    -2: 0.0010906577230485718,
    -1: 0.0022049272225017286,
    1: 0.020051104853608722,
    2: 0.0954139944889494,
    3: 0.004116641378803191,
    4: 0.002180124185625795,
    5: 6.922722705177504e-05,
    21: 0.007604124516892057}
"""
import logging
import numpy as np
import lhapdf

log = logging.getLogger(__name__)


class LHAPDFSet:
    """Wrapper for the lhapdf python interface.

    Once instantiated this class will load the PDF set from LHAPDF.
    If it is a T0 set only the CV will be loaded.

    For Monte Carlo sets the central value (member 0) is by default not included when taking
    the results for all members (i.e., when using ``grid_values``).
    """

    def __init__(self, name, error_type):
        log.info("PDF: %s ErrorType: %s", name, error_type)
        self._name = name
        self._error_type = error_type
        if self.is_t0:
            # If at this point we already know this is a T0 set, load only the CV
            self._lhapdf_set = [lhapdf.mkPDF(name)]
        else:
            self._lhapdf_set = lhapdf.mkPDFs(name)
        self._flavors = None

    @property
    def is_monte_carlo(self):
        """Check whether the error type is MC"""
        return self._error_type == "replicas"

    @property
    def is_t0(self):
        """Check whether we are in t0 mode"""
        return self._error_type == "t0"

    @property
    def n_members(self):
        """Return the number of active members in the PDF set"""
        return len(self.members)

    @property
    def members(self):
        """Return the members of the set, this depends on the error type:
        t0: returns only member 0
        MC: skip member 0
        Hessian: return all
        """
        if self.is_t0:
            return self._lhapdf_set[0:1]
        if self.is_monte_carlo:
            return self._lhapdf_set[1:]
        return self._lhapdf_set

    @property
    def central_member(self):
        """Returns a reference to member 0 of the PDF list"""
        return self._lhapdf_set[0]

    def xfxQ(self, x, q, member_idx, flavour):
        """Return the PDF value for one single point for one single member"""
        member_pdf = self.members[member_idx]
        res = member_pdf.xfxQ(x, q)
        return res[flavour]

    @property
    def flavors(self):
        """Returns the list of accepted flavors by the LHAPDF set"""
        if self._flavors is None:
            self._flavors = self.members[0].flavors()
        return self._flavors

    def grid_values(self, flavors: np.ndarray, xgrid: np.ndarray, qgrid: np.ndarray):
        """Returns the PDF values for every member for the required
        flavours, points in x and pointx in q
        The return shape is
            (members, flavors, xgrid, qgrid)
        Return
        ------
            ndarray of shape (members, flavors, xgrid, qgrid)
        Examples
        --------
        >>> import numpy as np
        >>> from validphys.lhapdfset import LHAPDFSet, ER_MC
        >>> pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", ER_MC)
        >>> xgrid = np.random.rand(10)
        >>> qgrid = np.random.rand(3)
        >>> flavs = np.arange(-4,4)
        >>> flavs[4] = 21
        >>> results = pdf.grid_values(flavs, xgrid, qgrid)
        """
        # Grid values loop
        ret = np.array(
            [
                [[member.xfxQ(flavors, x, q) for x in xgrid] for q in qgrid]
                for member in self.members
            ]
        )

        # Finally return in the documented shape
        return np.transpose(ret, (0, 3, 2, 1))
