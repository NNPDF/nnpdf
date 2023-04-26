"""
    Module containing an LHAPDF class compatible with validphys
    using the official lhapdf python interface.

    The ``.members`` and ``.central_member`` of the ``LHAPDFSet`` are
    LHAPDF objects (the typical output from ``mkPDFs``) and can be used normally.

    Examples
    --------
    >>> from validphys.lhapdfset import LHAPDFSet
    >>> pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", "replicas")
    >>> len(pdf.members)
    101
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
    """

    def __init__(self, name, error_type):
        self._name = name
        self._error_type = error_type
        if self.is_t0:
            # If at this point we already know this is a T0 set, load only the CV
            self._lhapdf_set = [lhapdf.mkPDF(name)]
        else:
            self._lhapdf_set = lhapdf.mkPDFs(name)
        self._flavors = None

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
        """Return the members of the set
        the special error type t0 returns only member 0
        """
        if self.is_t0:
            return self._lhapdf_set[0:1]
        return self._lhapdf_set

    @property
    def central_member(self):
        """Returns a reference to member 0 of the PDF list"""
        return self._lhapdf_set[0]

    def xfxQ(self, x, Q, n, fl):
        """Return the PDF value for one single point for one single member
        If the flavour is not included in the PDF (for instance top/antitop) return 0.0
        """
        member_pdf = self.members[n]
        res = member_pdf.xfxQ(x, Q)
        return res.get(fl, 0.0)

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
        >>> from validphys.lhapdfset import LHAPDFSet
        >>> pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", "replicas")
        >>> xgrid = np.random.rand(10)
        >>> qgrid = np.random.rand(3)
        >>> flavs = np.arange(-4,4)
        >>> flavs[4] = 21
        >>> results = pdf.grid_values(flavs, xgrid, qgrid)
        """
        # Create an array of x and q of equal length for LHAPDF
        xarr, qarr = (g.ravel() for g in np.meshgrid(xgrid, qgrid))
        # Ask LHAPDF for the values and swap the flavours and xgrid-qgrid axes
        raw = np.array([member.xfxQ(flavors, xarr, qarr) for member in self.members]).swapaxes(1, 2)
        # Unroll the xgrid-qgrid axes
        return raw.reshape(self.n_members, len(flavors), len(xgrid), len(qgrid))
