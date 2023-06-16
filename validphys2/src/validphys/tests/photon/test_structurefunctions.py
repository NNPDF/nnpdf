import numpy as np
import pineappl

from validphys.api import API
from validphys.core import PDF as PDFset
from validphys.lhapdfset import LHAPDFSet
import validphys.photon.structure_functions as sf

from ..conftest import PDF

TEST_THEORY = API.theoryid(theoryid=398)
NNPDF40 = PDFset(PDF)


class ZeroPdfs:
    def xfxQ(self, x, Q):
        res = {}
        for i in range(1, 6 + 1):
            res[i] = res[-i] = 0.0
        return res


class ZeroPdfs:
    def xfxQ(self, x, Q):
        res = {}
        for i in range(1, 6 + 1):
            res[i] = res[-i] = 0.0
        return res


def test_zero_pdfs():
    pdfs = ZeroPdfs()

    fake_theory = {
        "mc": 1.3,
        "mb": 5.0,
        "mt": 172.0,
        "kcThr": 1.0,
        "kbThr": 1.0,
        "ktThr": 1.0,
        "MaxNfPdf": 5,
    }

    f2lo = sf.F2LO(pdfs, fake_theory)

    np.testing.assert_equal(f2lo.thresh_t, np.inf)

    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(f2lo.fxq(x, Q), 0.0)


class FakeSet:
    def get_entry(self, string):
        return 0


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


class OnePdf:
    def __init__(self):
        self.ao = 1

    def xfxQ(self, x, Q):
        return 1.0

    def xfxQ2(self, x, Q):
        return 1.0**2

    def set(self):
        return FakeSet()


class OneFKTable:
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


class ZeroPdf:
    def __init__(self):
        self.ao = 1

    def xfxQ(self, x, Q):
        return 1.0

    def xfxQ2(self, x, Q):
        return 1.0**2

    def set(self):
        return FakeSet()


def test_F2(monkeypatch):
    # grid put to 0 and pdf put to 1
    monkeypatch.setattr(pineappl.fk_table.FkTable, "read", ZeroFKTable)
    structurefunc = sf.InterpStructureFunction("", OnePdf())
    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(structurefunc.fxq(x, Q), 0.0, rtol=1e-5)

    # grid put to 0 and real pdf
    pdf = LHAPDFSet("NNPDF40_nnlo_as_01180", "replicas")
    structurefunc = sf.InterpStructureFunction("", pdf.central_member)
    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(structurefunc.fxq(x, Q), 0.0, rtol=1e-5)

    # grid put to 1 and pdf put to 0
    monkeypatch.setattr(pineappl.fk_table.FkTable, "read", OneFKTable)
    structurefunc = sf.InterpStructureFunction("", ZeroPdf())
    for x in np.geomspace(1e-4, 1.0, 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(structurefunc.fxq(x, Q), 0.0, rtol=1e-5)


def test_params():
    pdfs = NNPDF40.load()
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
    pdfs = NNPDF40.load()
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
