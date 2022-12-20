from pathlib import Path
import pineappl
import numpy as np
from scipy.interpolate import RectBivariateSpline

class StructureFunction :
    def __init__(self, path_to_fktable, pdfs):
        self.path_to_fktable = Path(path_to_fktable)
        self.fktable = pineappl.fk_table.FkTable.read(path_to_fktable)
        self.pdfs = pdfs
        self.pdgid = int(pdfs.set().get_entry("Particle"))
        self.produce_interpolator()
    

    def produce_interpolator(self):
        x = np.unique(self.fktable.bin_left(1))
        q2 = np.unique(self.fktable.bin_left(0))
        predictions = self.fktable.convolute_with_one(self.pdgid, self.pdfs.xfxQ2)
        # here we require that the (x,Q2) couples that we passed
        # to pinefarm is a rectangular matrix

        grid2D = predictions.reshape(len(x),len(q2))
        # RectBivariateSpline is faster than interp2d
        self.interpolator = RectBivariateSpline(x, q2, grid2D)
    
    def FxQ(self, x, Q):
        return self.interpolator(x, Q**2)[0,0]

class F2LO :
    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
        self.Qmc = theory["Qmc"]
        self.Qmb = theory["Qmb"]
        self.Qmt = theory["Qmt"]
        if theory["MaxNfPdf"] <= 5 :
            self.Qmt = np.inf
        if theory["MaxNfPdf"] <= 4 :
            self.Qmb = np.inf
        if theory["MaxNfPdf"] <= 3 :
            self.Qmc = np.inf
        eu2 = 4. / 9
        ed2 = 1. / 9
        self.eq2 = [ed2, eu2, ed2, eu2, ed2, eu2] # d u s c b t
    
    def FxQ(self, x, Q):
        r"""
        Compute the LO DIS structure function F2.

        Parameters
        ----------
        x : float
            Bjorken's variable
        Q : float
            DIS hard scale
        
        Returns
        -------
        F2_LO : float
            Structure function F2 at LO
        """
        # at LO we use ZM-VFS
        if Q < self.Qmc :
            nf = 3
        elif Q < self.Qmb :
            nf = 4
        elif Q < self.Qmt :
            nf = 5
        else :
            nf = 6
        res = 0
        for i in range(1, nf+1):
            res += self.eq2[i-1] * (self.pdfs.xfxQ(x, Q)[i] + self.pdfs.xfxQ(x, Q)[-i])
        return res
