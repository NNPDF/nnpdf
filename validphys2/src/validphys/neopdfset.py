"""
Module containing the NeoPDF-backend PDF set compatible with validphys.

``NeoPDFSet`` has the same public interface as ``LHAPDFSet`` so it can
be used as a drop-in replacement anywhere an ``LHAPDFSet`` is expected.

Examples
--------
>>> from validphys.neopdfset import NeoPDFSet
>>> pdf = NeoPDFSet("NNPDF40_nnlo_as_01180", "replicas")
>>> len(pdf.members)
101
>>> pdf.grid_values([-2, -1, 1, 2, 21], [0.1, 0.5], [10.0, 100.0]).shape
(101, 5, 2, 2)
"""

import logging

from neopdf.pdf import PDF as _NeoPDF
import numpy as np

log = logging.getLogger(__name__)


class NeoPDFSet:
    """Wrapper around the NeoPDF interpolation library.

    Provides the same interface as ``LHAPDFSet`` so it can be used
    interchangeably within validphys and n3fit.
    """

    def __init__(self, name: str, error_type: str):
        self._name = name
        self._error_type = error_type
        self._flavors = None
        if self.is_t0:
            self._members = [_NeoPDF(name)]
        else:
            self._members = _NeoPDF.mkPDFs(name)

    @property
    def is_t0(self) -> bool:
        """Check whether it is a t0 set."""
        return self._error_type == "t0"

    @property
    def n_members(self) -> int:
        """Return the total members of the PDF set."""
        return len(self.members)

    @property
    def members(self):
        """Return a given member of the PDF set.

        The special error type t0 returns only member 0.
        """
        if self.is_t0:
            return self._members[0:1]
        return self._members

    @property
    def central_member(self):
        """Return a reference to member 0."""
        return self._members[0]

    def xfxQ(self, x: float, Q: float, n: int, fl: int) -> float:
        """Return the PDF value for a single point for a single member.

        If the flavour is absent from the set, simply return 0.0.
        """
        if fl not in self.flavors:
            return 0.0
        return self.members[n].xfxQ2(fl, x, Q**2)

    @property
    def flavors(self):
        """Return the list of particle IDs supported by this PDF set."""
        if self._flavors is None:
            self._flavors = self.members[0].pids()
        return self._flavors

    def grid_values(self, flavors: np.ndarray, xgrid: np.ndarray, qgrid: np.ndarray) -> np.ndarray:
        """Return PDF values on a grid for all the members.

        Parameters
        ----------
        flavors:
            Array of PDG particle IDs.
        xgrid:
            Array of x values.
        qgrid:
            Array of Q values (not Q²).

        Returns
        -------
        ndarray of shape ``(members, flavors, len(xgrid), len(qgrid))``
        """
        q2grid = np.asarray(qgrid) ** 2
        xgrid = np.asarray(xgrid)
        pids = list(flavors)
        nx, nq, nfl = len(xgrid), len(qgrid), len(pids)

        results = []
        for member in self.members:
            vals = np.array(member.xfxQ2s(pids, xgrid, q2grid))
            results.append(vals.reshape(nfl, nx, nq))
        return np.array(results)
