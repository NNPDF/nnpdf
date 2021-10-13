"""
    Module containing an LHAPDF class compatible with validphys
    using the official lhapdf python interface.

    It exposes an interface (almost) compatible with libNNPDF::LHAPDFSet
"""
import logging
from typing import NamedTuple
import numpy as np
import lhapdf

log = logging.getLogger(__name__)


class PDFErrorType(NamedTuple):
    """Namedtuple containing the information about the error type of the PDF"""

    name: str
    monte_carlo: bool
    libNNPDF: int
    t0: bool
    description: str


ER_NONE = PDFErrorType("erType_ER_NONE", False, 0, False, "No error info")
ER_MC = PDFErrorType("erType_ER_MC", True, 1, False, "Monte Carlo")
ER_MC68 = PDFErrorType("erType_ER_MC68", True, 2, False, "Monte Carlo 68pc")
ER_MCT0 = PDFErrorType("erType_ER_MCT0", False, 3, True, "Monte Carlo t0")
ER_EIG = PDFErrorType("erType_ER_EIG", False, 4, False, "Eigenvector 68pc")
ER_EIG90 = PDFErrorType("erType_ER_EIG90", False, 5, False, "Eigenvector 90pc")
ER_SYMEIG = PDFErrorType("erType_ER_SYMEIG", False, 6, False, "Symmetric eigenvectors")
_libNNPDF_errors = [ER_NONE, ER_MC, ER_MC68, ER_MCT0, ER_EIG, ER_EIG90, ER_SYMEIG]


def _parse_flavs(f):
    # TODO: for now we believe what the user unless they say 0, then they're wrong
    if f == 0:
        return 21
    return f


class LHAPDFSet:
    """Wrapper for the lhapdf python interface"""

    def __init__(self, name, error_type, lhapdf_verbosity=0):
        if isinstance(error_type, int):
            # libNNPDF error types were int
            error_type = _libNNPDF_errors[error_type]
        self._name = name
        self._error_type = error_type
        self._flavors = None
        self._lhapdf_set = lhapdf.mkPDFs(name)
        self._libNNPDF_set = None
        self.legacy_interface()
        log.info("PDF: %s ErrorType: %s", name, error_type.name)
        lhapdf.setVerbosity(lhapdf_verbosity)

    def legacy_interface(self):
        """Setting some methods and attribute as per libNNPDF specs"""
        self.GetMembers = self.get_members
        self.GetSetName = lambda self: self.name

    def get_members(self):
        """Get the number of members"""
        return len(self.members)

    @property
    def is_monte_carlo(self):
        """Check whether the error type is MC"""
        return self._error_type.monte_carlo

    @property
    def is_t0(self):
        """Check whether we are in t0 mode"""
        return self._error_type.t0

    @property
    def members(self):
        """Return the members of the set, this depends on the error type:
        t0: returns only member 0
        MC: skip member 0
        Eigen: return all
        """
        if self.is_t0:
            return self._lhapdf_set[0:1]
        if self.is_monte_carlo:
            return self._lhapdf_set[1:]
        return self._lhapdf_set

    def as_libNNPDF(self):
        """Imports libNNPDF to return the PDF set in NNPDF form
        needed for classes that still rely in libNNPDF and will need a PDF
        (at the moment due to closure tests)
        """
        if self._libNNPDF_set is None:
            # TODO obvs
            from NNPDF import LHAPDFSet

            self._libNNPDF_set = LHAPDFSet(self._name, self._error_type.libNNPDF)
        return self._libNNPDF_set

    def xfxQ(self, x, q, member_idx, flavour):
        """Return the PDF value for one single point"""
        member_pdf = self.members[member_idx]
        res = member_pdf.xfxQ(x, q)
        return res[_parse_flavs(flavour)]

    @property
    def flavors(self):
        """Returns the list of accepted flavors by the LHAPDF set"""
        if self._flavors is None:
            self._flavors = self.members[0].flavors()
        return self._flavors

    def grid_values(self, flavors, xgrid, qgrid):
        """Reimplementation of libNNPDF grid_values
        The return shape is
            (members, flavors, xgrid, qgrid)

        Examples
        --------
        >>> import numpy as np
        >>> from validphys.lhapdfset import LHAPDFSet, ER_MC
        >>> pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", ER_MC)
        >>> xgrid = np.random.rand(10)
        >>> qgrid = np.random.rand(3)
        >>> flavs = np.arange(-4,4)
        >>> results = valid_pdf.grid_values(flavs, xgrid, qgrid)
        """
        # TODO Digest the input flavors in a less confusing way
        parsed_flavors = [_parse_flavs(f) for f in flavors]
        flavs = []
        non_flavs = []
        for pf in parsed_flavors:
            if pf in self.flavors:
                flavs.append(pf)
            else:
                non_flavs.append(pf)
        xvals = []
        for x in xgrid:
            qvals = []
            for q in qgrid:
                mres = []
                for member in self.members:
                    all_flavs_dict = member.xfxQ(x, q)
                    mres.append([all_flavs_dict[i] for i in flavs])
                qvals.append(mres)
            xvals.append(qvals)
        ret = np.array(xvals)
        # Add 0s for flavours that were not in the PDF but the user asked for
        # Could equally well be done in the loop above but looked like an unnecesary complication
        zeros = np.zeros_like(ret[..., 0])
        for f in non_flavs:
            idx = parsed_flavors.index(f)
            # Do it one by one to keep the index in the right places
            ret = np.insert(ret, idx, zeros, -1)
        return np.transpose(ret, (2, 3, 0, 1))


for error_type in _libNNPDF_errors:
    setattr(LHAPDFSet, error_type.name, error_type)
