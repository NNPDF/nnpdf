import pytest
from validphys.photon.compute import Photon
from validphys.photon import structure_functions
import lhapdf
import fiatlux
import numpy as np
from collections import namedtuple

class faketheory():
    def get_description(self):
        return {
            "alphaqed": 0.01,
            "Qref": 91.2,
            "Qmc": 1.3,
            "Qmb": 4.92,
            "Qmt": 173.,
            "MaxNfAs": 5,

        }

fiatlux_runcard = {
    "pdf_name": "no_pdf",
    "path_to_F2": "/",
    "path_to_FL": "/",
    "additional_errors": False, 
}

photon = namedtuple('photon', ['total', 'elastic', 'inelastic'])

class fakeFiatlux():
    def __init__(self, runcard):
        self.runcard = runcard
        self.alphaem = None
        self.qref = None
        self.trash1 = None
        self.trash2 = None
        self.f2 = None
        self.fl = None
        self.f2lo = None
        self.res = photon(0, 0, 0)
    
    def PlugAlphaQED(self, alphaem, qref):
        self.alphaem = alphaem
        self.qref = qref
    
    def InsertInelasticSplitQ(self, args):
        self.trash1 = args[0]
        self.trash2 = args[1]
    
    def PlugStructureFunctions(self, f2, fl, f2lo):
        self.f2 = f2
        self.fl = fl
        self.f2lo = f2lo
    
    def EvaluatePhoton(self, x, q):
        return self.res
    
class fakeStructureFunction():
    def __init__(self, path, pdfs):
        self.path = path
        self.pdfs = pdfs
    
    def FxQ(self):
        return 0

class fakeF2LO():
    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
        self.theory = theory
    
    def FxQ(self):
        return 0


def test_init(monkeypatch):
    monkeypatch.setattr(structure_functions, "StructureFunction", fakeStructureFunction)
    monkeypatch.setattr(structure_functions, "F2LO", fakeF2LO)
    monkeypatch.setattr(lhapdf, "mkPDF", lambda *args: 1)
    monkeypatch.setattr(fiatlux, "FiatLux", fakeFiatlux)
    monkeypatch.setattr(Photon, "produce_interpolators", lambda *args: None)
    photon = Photon(faketheory(), fiatlux_runcard, [1,2,3])
    np.testing.assert_equal(photon.replicas_id, [1,2,3])
    np.testing.assert_equal(photon.Qmt, np.inf)
    np.testing.assert_almost_equal(photon.Qmb, 4.92)
    np.testing.assert_almost_equal(photon.Qmc, 1.3)
    np.testing.assert_almost_equal(photon.thresh[5], 91.2)
    np.testing.assert_almost_equal(photon.thresh[4], 4.92)
    np.testing.assert_almost_equal(photon.thresh[3], 1.3)
    np.testing.assert_almost_equal(photon.alpha_thresh[5], 0.01)
    np.testing.assert_almost_equal(photon.alpha_thresh[4], photon.alpha_em_nlo(4.92, 0.01, 91.2, 5))
    np.testing.assert_almost_equal(photon.alpha_thresh[3], photon.alpha_em_nlo(1.3, photon.alpha_thresh[4], 4.92, 4))
    np.testing.assert_equal(len(photon.alpha_thresh), 3)
    np.testing.assert_equal(len(photon.thresh), 3)


