from collections import namedtuple
from pathlib import Path

import fiatlux
import numpy as np
from validphys.photon import structure_functions
from validphys.photon.compute import Photon

from ..conftest import PDF


class FakeTheory:
    def __init__(self):
        self.path = Path("/fake/path/")

    def get_description(self):
        return {
            "alphaqed": 0.01,
            "Qref": 91.2,
            "mc": 1.3,
            "mb": 4.92,
            "mt": 172.0,
            "kcThr": 1.0,
            "kbThr": 1.0,
            "ktThr": 1.0,
            "MaxNfAs": 5,
            "MaxNfPdf": 5,
        }


fiatlux_runcard = {
    "pdf_name": PDF,
    "additional_errors": False,
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
    monkeypatch.setattr(Photon, "produce_interpolators", lambda *args: None)

    photon = Photon(FakeTheory(), fiatlux_runcard, [1, 2, 3])

    np.testing.assert_equal(photon.replicas_id, [1, 2, 3])
    np.testing.assert_equal(photon.fiatlux_runcard, fiatlux_runcard)
    np.testing.assert_almost_equal(photon.q_in2, 1e4)
    np.testing.assert_almost_equal(
        photon.alpha_em_ref, FakeTheory().get_description()["alphaqed"]
    )

def test_masses_init(monkeypatch):
    monkeypatch.setattr(
        structure_functions, "InterpStructureFunction", FakeStructureFunction
    )
    monkeypatch.setattr(structure_functions, "F2LO", FakeF2LO)

    monkeypatch.setattr(fiatlux, "FiatLux", FakeFiatlux)
    monkeypatch.setattr(Photon, "produce_interpolators", lambda *args: None)
    photon = Photon(FakeTheory(), fiatlux_runcard, [1, 2, 3])
    np.testing.assert_equal(photon.thresh_t, np.inf)
    np.testing.assert_almost_equal(photon.thresh_b, 4.92)
    np.testing.assert_almost_equal(photon.thresh_c, 1.3)

def test_set_thresholds_alpha_em(monkeypatch):
    monkeypatch.setattr(
        structure_functions, "InterpStructureFunction", FakeStructureFunction
    )
    monkeypatch.setattr(structure_functions, "F2LO", FakeF2LO)

    monkeypatch.setattr(fiatlux, "FiatLux", FakeFiatlux)
    monkeypatch.setattr(Photon, "produce_interpolators", lambda *args: None)
    photon = Photon(FakeTheory(), fiatlux_runcard, [1, 2, 3])
    np.testing.assert_almost_equal(photon.thresh[5], 91.2)
    np.testing.assert_almost_equal(photon.thresh[4], 4.92)
    np.testing.assert_almost_equal(photon.thresh[3], 1.3)
    np.testing.assert_almost_equal(photon.alpha_thresh[5], 0.01)
    np.testing.assert_almost_equal(
        photon.alpha_thresh[4], photon.alpha_em_nlo(4.92, 0.01, 91.2, 5)
    )
    np.testing.assert_almost_equal(
        photon.alpha_thresh[3],
        photon.alpha_em_nlo(1.3, photon.alpha_thresh[4], 4.92, 4),
    )
    np.testing.assert_equal(len(photon.alpha_thresh), 3)
    np.testing.assert_equal(len(photon.thresh), 3)

def test_betas(monkeypatch):
    monkeypatch.setattr(
        structure_functions, "InterpStructureFunction", FakeStructureFunction
    )
    monkeypatch.setattr(structure_functions, "F2LO", FakeF2LO)

    monkeypatch.setattr(fiatlux, "FiatLux", FakeFiatlux)
    monkeypatch.setattr(Photon, "produce_interpolators", lambda *args: None)
    photon = Photon(FakeTheory(), fiatlux_runcard, [1, 2, 3])
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
        np.testing.assert_allclose(photon.beta0[nf], vec_beta0[nf - 3], rtol=1e-7)
        np.testing.assert_allclose(photon.b1[nf], vec_b1[nf - 3], rtol=1e-7)
