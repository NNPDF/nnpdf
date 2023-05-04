from collections import namedtuple
from pathlib import Path

import fiatlux
import numpy as np
from validphys.photon import structure_functions
from validphys.photon.compute import Photon, Alpha
from validphys.core import PDF as PDFset

from ..conftest import PDF


class FakeTheory:
    def __init__(self):
        self.path = Path("/fake/path/")

    def get_description(self):
        return {
            "alphaqed": 0.01,
            "Qref": 91.2,
            "QrefQED": 91.2,
            "mc": 1.3,
            "mb": 4.92,
            "mt": 172.0,
            "kcThr": 1.0,
            "kbThr": 1.0,
            "ktThr": 1.0,
            "MaxNfAs": 5,
            "MaxNfPdf": 5,
            "MP": 0.938
        }


fiatlux_runcard = {
    "luxset": PDFset(PDF),
    "additional_errors": False,
    "luxseed": 123456789
}

photon = namedtuple("photon", ["total", "elastic", "inelastic"])


class FakeFiatlux:
    def __init__(self, runcard):
        self.runcard = runcard
        self.alphaem = None
        self.qref = None
        self.trash1 = None
        self.trash2 = None
        self.f2 = None
        self.fl = None
        self.f2lo = None
        self.res = photon(0., 0., 0.)

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


class FakeStructureFunction:
    def __init__(self, path, pdfs):
        self.path = path
        self.pdfs = pdfs

    def fxq(self):
        return 0


class FakeF2LO:
    def __init__(self, pdfs, theory):
        self.pdfs = pdfs
        self.theory = theory

    def fxq(self):
        return 0


def test_parameters_init(monkeypatch):
    monkeypatch.setattr(
        structure_functions, "InterpStructureFunction", FakeStructureFunction
    )
    monkeypatch.setattr(structure_functions, "F2LO", FakeF2LO)
    monkeypatch.setattr(fiatlux, "FiatLux", FakeFiatlux)
    monkeypatch.setattr(Photon, "compute_photon_array", lambda *args: np.zeros(196))

    photon = Photon(FakeTheory(), fiatlux_runcard, [1, 2, 3])
    alpha = Alpha(FakeTheory().get_description())

    np.testing.assert_equal(photon.replicas, [1, 2, 3])
    np.testing.assert_equal(photon.luxpdfset._name, fiatlux_runcard["luxset"].name)
    np.testing.assert_equal(photon.additional_errors, fiatlux_runcard["additional_errors"])
    np.testing.assert_equal(photon.luxseed, fiatlux_runcard["luxseed"])
    np.testing.assert_almost_equal(
        alpha.alpha_em_ref, FakeTheory().get_description()["alphaqed"]
    )

def test_masses_init():
    alpha = Alpha(FakeTheory().get_description())
    np.testing.assert_equal(alpha.thresh_t, np.inf)
    np.testing.assert_almost_equal(alpha.thresh_b, 4.92)
    np.testing.assert_almost_equal(alpha.thresh_c, 1.3)

def test_set_thresholds_alpha_em(monkeypatch):
    monkeypatch.setattr(
        structure_functions, "InterpStructureFunction", FakeStructureFunction
    )
    monkeypatch.setattr(structure_functions, "F2LO", FakeF2LO)

    monkeypatch.setattr(fiatlux, "FiatLux", FakeFiatlux)
    monkeypatch.setattr(Photon, "compute_photon_array", lambda *args: np.zeros(196))

    alpha = Alpha(FakeTheory().get_description())

    np.testing.assert_almost_equal(alpha.thresh[5], 91.2)
    np.testing.assert_almost_equal(alpha.thresh[4], 4.92)
    np.testing.assert_almost_equal(alpha.thresh[3], 1.3)
    np.testing.assert_almost_equal(alpha.alpha_thresh[5], 0.01)
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[4], alpha.alpha_em_fixed_flavor(4.92, 0.01, 91.2, 5)
    )
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[3],
        alpha.alpha_em_fixed_flavor(1.3, alpha.alpha_thresh[4], 4.92, 4),
    )
    np.testing.assert_equal(len(alpha.alpha_thresh), 3)
    np.testing.assert_equal(len(alpha.thresh), 3)

def test_betas():
    alpha = Alpha(FakeTheory().get_description())
    vec_beta0 = [
        -0.5305164769729844,
        -0.6719875374991137,
        -0.7073553026306458,
        -0.8488263631567751,
    ]
    vec_b1 = [
        0.17507043740108488,
        0.1605510390839295,
        0.1538497783221655,
        0.1458920311675707,
    ]
    for nf in range(3, 6 + 1):
        np.testing.assert_allclose(alpha.beta0[nf], vec_beta0[nf - 3], rtol=1e-7)
        np.testing.assert_allclose(alpha.b1[nf], vec_b1[nf - 3], rtol=1e-7)
