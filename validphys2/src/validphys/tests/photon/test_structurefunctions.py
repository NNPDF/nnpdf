import numpy as np
import pineappl

from validphys.api import API
from validphys.core import PDF as PDFset
import validphys.photon.structure_functions as sf

from ..conftest import PDF

TEST_THEORY = API.theoryid(theoryid=398)


class ZeroPdfs:
    def xfxQ(self, x, Q):
        res = {}
        for i in range(1, 6 + 1):
            res[i] = res[-i] = 0.0
        res[21] = 0.0
        res[22] = 0.0
        return res
    
    def xfxQ2(self, i, x, Q2):
        return self.xfxQ(x, np.sqrt(Q2))[i]


class ZeroFKTable:
    def __init__(self, path):
        self.path = path
        self.xgrid = np.geomspace(1e-4, 1.0, 10)
        self.qgrid = np.geomspace(1.65, 1000, 10)

    def bin_left(self, i):
        if i == 1:
            return self.xgrid
        if i == 0:
            return self.qgrid
        else:
            return 0

    def convolute_with_one(self, pdgid, xfxQ2):
        return np.zeros((10, 10))


def test_zero_pdfs():
    "test that a zero PDF gives a zero structure function"
    pdfs = ZeroPdfs()
    theory = TEST_THEORY.get_description()
    path_to_F2 = TEST_THEORY.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
    path_to_FL = TEST_THEORY.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"

    f2 = sf.InterpStructureFunction(path_to_F2, pdfs)
    fl = sf.InterpStructureFunction(path_to_FL, pdfs)
    f2lo = sf.F2LO(pdfs, theory)

    np.testing.assert_equal(f2lo.thresh_t, np.inf)

    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(f2lo.fxq(x, Q), 0.0)
            np.testing.assert_allclose(f2.fxq(x, Q), 0.0)
            np.testing.assert_allclose(fl.fxq(x, Q), 0.0)



def test_zero_grid(monkeypatch):
    "test that a zero grid gives a zero structure function"
    # patching pineappl.fk_table.FkTable to use ZeroFKTable
    monkeypatch.setattr(pineappl.fk_table.FkTable, "read", ZeroFKTable)
    pdfs = PDFset(PDF).load()
    structurefunc = sf.InterpStructureFunction("", pdfs.central_member)
    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(structurefunc.fxq(x, Q), 0.0, rtol=1e-5)



def test_params():
    "test initialization of parameters"
    pdfs = PDFset(PDF).load()
    replica = 1
    theory = TEST_THEORY.get_description()
    for channel in ["F2", "FL"]:
        tmp = "fastkernel/fiatlux_dis_" + channel + ".pineappl.lz4"
        path_to_fktable = TEST_THEORY.path / tmp
        struct_func = sf.InterpStructureFunction(path_to_fktable, pdfs.members[replica])
        np.testing.assert_allclose(struct_func.q2_max, 1e8)
    f2lo = sf.F2LO(pdfs.members[replica], theory)
    np.testing.assert_allclose(f2lo.thresh_c, theory["mc"])
    np.testing.assert_allclose(f2lo.thresh_b, theory["mb"])
    np.testing.assert_allclose(f2lo.thresh_t, np.inf)


def test_interpolation_grid():
    """test that the values coming out of InterpStructureFunction match the grid ones"""
    pdfs = PDFset(PDF).load()
    for replica in [1, 2, 3]:
        for channel in ["F2", "FL"]:
            tmp = "fastkernel/fiatlux_dis_" + channel + ".pineappl.lz4"
            path_to_fktable = TEST_THEORY.path / tmp
            fktable = pineappl.fk_table.FkTable.read(path_to_fktable)
            x = np.unique(fktable.bin_left(1))
            q2 = np.unique(fktable.bin_left(0))
            predictions = fktable.convolute_with_one(2212, pdfs.members[replica].xfxQ2)
            grid2D = predictions.reshape(len(x), len(q2))

            struct_func = sf.InterpStructureFunction(path_to_fktable, pdfs.members[replica])
            for i, x_ in enumerate(x):
                for j, q2_ in enumerate(q2):
                    np.testing.assert_allclose(
                        struct_func.fxq(x_, np.sqrt(q2_)), grid2D[i, j], rtol=1e-5
                    )
