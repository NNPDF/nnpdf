import pineappl
import numpy as np
from scipy.interpolate import RectBivariateSpline

class StructureFunction :
    def __init__(self, path_to_fktable, pdfs):
        self.fktable = pineappl.fk_table.FkTable.read(path_to_fktable)
        self.pdfs = pdfs
        self.pdgid = int(pdfs.set().get_entry("Particle"))
        self.produce_interpolator()
    

    def produce_interpolator(self):
        x = np.unique(self.fktable.bin_left(1))
        q2 = np.unique(self.fktable.bin_left(0))
        self.xmin = min(x)
        self.qmin = min(np.sqrt(q2))
        self.qmax = max(np.sqrt(q2))

        predictions = self.fktable.convolute_with_one(self.pdgid, self.pdfs.xfxQ2)
        # here we require that the (x,Q2) couples that we passed
        # to pinefarm is a rectangular matrix

        grid2D = predictions.reshape(len(x), len(q2))
        self.interpolator = RectBivariateSpline(x, q2, grid2D)
    
    def FxQ(self, x, Q):
        # here we are requiring that the grid that we pass to fiatlux
        # has Qmin = 1 (fiatlux doesn't go below Q=1)
        if x < self.xmin or Q > self.qmax :
            return 0.
        return self.interpolator(x, Q**2)[0, 0]

class F2LO :
    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
        # TODO : maybe they shoud be kDIS instead of k, but usually they are the same
        self.thresh_c = theory["kcThr"] * theory["mc"]
        self.thresh_b = theory["kbThr"] * theory["mb"]
        self.thresh_t = theory["ktThr"] * theory["mt"]
        if theory["MaxNfPdf"] <= 5 :
            self.thresh_t = np.inf
        if theory["MaxNfPdf"] <= 4 :
            self.thresh_b = np.inf
        if theory["MaxNfPdf"] <= 3 :
            self.thresh_c = np.inf
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
        if Q < self.thresh_c :
            nf = 3
        elif Q < self.thresh_b :
            nf = 4
        elif Q < self.thresh_t :
            nf = 5
        else :
            nf = 6
        res = 0
        pdfs_values = self.pdfs.xfxQ(x, Q)
        for i in range(1, nf+1):
            res += self.eq2[i-1] * (pdfs_values[i] + pdfs_values[-i])
        return res
