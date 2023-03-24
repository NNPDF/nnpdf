import validphys.photon.structure_functions as sf
import numpy as np
import pineappl

def test_zero_pdfs():
    class fake_pdfs:
        def xfxQ(self, x, Q):
            res = {}
            for i in range(1, 6+1):
                res[i] = res[-i] = 0.
            return res

    pdfs = fake_pdfs()

    fake_theory = {
        "mc": 1.3,
        "mb": 5. ,
        "mt": 172.,
        "kcThr": 1.,
        "kbThr": 1.,
        "ktThr": 1.,
        "MaxNfPdf": 5,
    }

    f2lo = sf.F2LO(pdfs, fake_theory)

    np.testing.assert_equal(f2lo.thresh_t, np.inf)

    for x in np.geomspace(1e-4, 1., 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(f2lo.fxq(x, Q), 0.)

class FakeFKTable():
    def __init__(self, path):
        self.path = path
        self.xgrid = np.geomspace(1e-4, 1., 10)
        self.qgrid = np.geomspace(1.65, 1000, 10)

    def bin_left(self, i):
        if i == 1:
            return self.xgrid
        if i == 0 :
            return self.qgrid
        else:
            return 0

    def convolute_with_one(self, pdgid, xfxQ2):
        return np.zeros((10, 10))

class fakeset():
    def get_entry(self, string):
        return 0

class fakepdf():
    def __init__(self):
        self.ao = 1
    
    def xfxQ(self, x, Q):
        return 1.
    
    def xfxQ2(self, x, Q):
        return 1.**2
    
    def set(self):
        return fakeset()


def test_F2(monkeypatch):
    monkeypatch.setattr(pineappl.fk_table.FkTable, "read", FakeFKTable)
    structurefunc = sf.StructureFunction("", fakepdf())
    for x in np.geomspace(1e-4, 1., 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(structurefunc.fxq(x, Q), 0., rtol=1e-5)

    