from abc import ABC, abstractmethod

import numpy as np
import pineappl
from scipy.interpolate import RectBivariateSpline

from . import constants


class StructureFunction(ABC):
    """Abstract class for the DIS structure functions"""

    @abstractmethod
    def fxq(self, x, q):
        r"""
        Abstract method returning the value of a DIS structure functions.

        Parameters
        ----------
        x : float
            Bjorken x
        q : float
            DIS hard scale

        Returns
        -------
        F_{2,L}: float
            Structure function F2 or FL
        """
        pass


class InterpStructureFunction(StructureFunction):
    """
    Compute an interpolated structure function convoluting an FKtable
    with a PDF.
    """

    def __init__(self, path_to_fktable, pdfs):
        self.fktable = pineappl.fk_table.FkTable.read(path_to_fktable)
        x = np.unique(self.fktable.bin_left(1))
        q2 = np.unique(self.fktable.bin_left(0))

        self.q2_max = max(q2)

        predictions = self.fktable.convolute_with_one(2212, pdfs.xfxQ2)

        grid2D = predictions.reshape(len(x), len(q2))
        self.interpolator = RectBivariateSpline(x, q2, grid2D)

    def fxq(self, x, q):
        r"""
        Compute the DIS structure function interpolating the grid.

        Parameters
        ----------
        x : float
            Bjorken x
        q : float
            DIS hard scale

        Returns
        -------
        F_{2,L}: float
            Structure function F2 or FL
        """
        return self.interpolator(x, q**2)[0, 0]


class F2LO(StructureFunction):
    """Analytically compute the structure function F2 at leading order."""

    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
        self.thresh_c = theory["kcThr"] * theory["mc"]
        self.thresh_b = theory["kbThr"] * theory["mb"]
        self.thresh_t = theory["ktThr"] * theory["mt"]
        if theory["MaxNfPdf"] <= 5:
            self.thresh_t = np.inf
        if theory["MaxNfPdf"] <= 4:
            self.thresh_b = np.inf
        if theory["MaxNfPdf"] <= 3:
            self.thresh_c = np.inf
        if theory["IC"] == 1:
            self.thresh_c = 0.0
        self.eq2 = [constants.ED2, constants.EU2] * 3  # d u s c b t

    def fxq(self, x, q):
        r"""
        Analytically compute F2LO.

        Parameters
        ----------
        x : float
            Bjorken x
        q : float
            DIS hard scale

        Returns
        -------
        F_2^{LO} : float
            Structure function F2 at LO
        """
        # at LO we use ZM-VFS
        if q < self.thresh_c:
            nf = 3
        elif q < self.thresh_b:
            nf = 4
        elif q < self.thresh_t:
            nf = 5
        else:
            nf = 6
        res = 0
        pdfs_values = self.pdfs.xfxQ(x, q)
        for i in range(1, nf + 1):
            res += self.eq2[i - 1] * (pdfs_values[i] + pdfs_values[-i])
        return res
