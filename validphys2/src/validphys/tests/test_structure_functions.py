import apfel
from validphys.photon.structure_functions import StructureFunction
import lhapdf
import numpy as np
from pathlib import Path
from eko.interpolation import make_grid
import pineappl
from scipy.interpolate import interp2d, RectBivariateSpline

pdf_name = "NNPDF40_nnlo_as_01180"
pdfs = lhapdf.mkPDF(pdf_name, 0)
pdgid = int(pdfs.set().get_entry("Particle"))
path_to_F2 = '/home/laurenti/n3pdf/pineline/data/fktables/200/DIS_F2.pineappl.lz4'
path_to_FL = '/home/laurenti/n3pdf/pineline/data/fktables/200/DIS_FL.pineappl.lz4'

fktablef2=pineappl.fk_table.FkTable.read(Path(path_to_F2))
fktablefl=pineappl.fk_table.FkTable.read(Path(path_to_FL))

def produce_interpolator1(fktable):
    x = np.unique(fktable.bin_left(1))
    q2 = np.unique(fktable.bin_left(0))
    predictions = fktable.convolute_with_one(pdgid, pdfs.xfxQ2)
    grid2D = predictions.reshape(len(x),len(q2)).T
    interpolator = interp2d(x, q2, grid2D)
    return interpolator

def produce_interpolator2(fktable):
    x = np.unique(fktable.bin_left(1))
    q2 = np.unique(fktable.bin_left(0))
    predictions = fktable.convolute_with_one(pdgid, pdfs.xfxQ2)
    grid2D = predictions.reshape(len(x),len(q2))
    print(grid2D)
    interpolator = RectBivariateSpline(x, q2, grid2D)
    return interpolator

inter1_f2 = produce_interpolator1(fktablef2)
inter2_f2 = produce_interpolator2(fktablef2)
inter1_fl = produce_interpolator1(fktablefl)
inter2_fl = produce_interpolator2(fktablefl)

for x in make_grid(10,10):
    for q in [5, 10, 20, 50, 100]:
        print(x, q)
        np.testing.assert_allclose(inter1_f2(x, q**2)[0], inter2_f2(x, q**2)[0,0], rtol=21e-1)
        np.testing.assert_allclose(inter1_fl(x, q**2)[0], inter2_fl(x, q**2)[0,0], rtol=21e-1)