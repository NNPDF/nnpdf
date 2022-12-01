from pathlib import Path
import pineappl
import numpy as np
from scipy.interpolate import interp2d
# from scipy.interpolate import RectBivariateSpline

class StructureFunction :
    def __init__(self, path_to_fktable, pdfs):
        # self.path_to_grid = Path(path_to_grid)
        self.path_to_fktable = Path(path_to_fktable)
        # self.grid = pineappl.grid.Grid.read(path_to_grid)
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
        grid2D = predictions.reshape(len(x),len(q2)).T
        # TODO: are len(x) and len(q2) in the correct order?
        self.interpolator = interp2d(x, q2, grid2D)
        # consider using scipy.RectBivariateSpline.
    
    def FxQ(self, x, Q):
        return self.interpolator(x, Q**2)


