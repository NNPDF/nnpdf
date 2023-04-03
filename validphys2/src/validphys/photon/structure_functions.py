import pineappl
import numpy as np
from scipy.interpolate import RectBivariateSpline

class StructureFunction :
    """Abstract class for the DIS structure functions"""
    def fxq(self, x, Q) :
        pass

class InterpStructureFunction(StructureFunction):
    """
    Compute an interpolated structure function convoluting an FKtable
    with a PDF.
    """
    def __init__(self, path_to_fktable, pdfs):
        self.fktable = pineappl.fk_table.FkTable.read(path_to_fktable)
        self.pdfs = pdfs
        self.pdgid = int(pdfs.set().get_entry("Particle"))
        self.produce_interpolator()
    

    def produce_interpolator(self):
        """Produce the interpolation function to be called in fxq."""
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
    
    def fxq(self, x, q):
        r"""
        Compute the DIS structure function interpolating the grid.

        Parameters
        ----------
        x : float
            Bjorken's variable
        Q : float
            DIS hard scale
        
        Returns
        -------
        F_{2,L}: float
            Structure function F2 or FL
        """
        # here we are requiring that the grid that we pass to fiatlux
        # has Qmin = 1 (fiatlux doesn't go below Q=1)
        # if x < self.xmin or q > self.qmax :
        #     return 0.
        return self.interpolator(x, q**2)[0, 0]

class F2LO(StructureFunction) :
    """Compute analytically the leading order structure function for F2."""
    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
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
    
    def fxq(self, x, q):
        r"""
        Compute the analytical form of F2LO.

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
        if q < self.thresh_c :
            nf = 3
        elif q < self.thresh_b :
            nf = 4
        elif q < self.thresh_t :
            nf = 5
        else :
            nf = 6
        res = 0
        pdfs_values = self.pdfs.xfxQ(x, q)
        for i in range(1, nf+1):
            res += self.eq2[i-1] * (pdfs_values[i] + pdfs_values[-i])
        return res
