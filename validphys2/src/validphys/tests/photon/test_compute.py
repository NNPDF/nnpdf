from collections import namedtuple
from pathlib import Path

import numpy as np
from validphys.photon import structure_functions
from validphys.photon.compute import Photon, Alpha
from validphys.core import PDF as PDFset
from validphys.api import API

from ..conftest import PDF
from eko.io import EKO

TEST_THEORY = API.theoryid(theoryid=398)

FIATLUX_RUNCARD = {
    "luxset": PDFset(PDF),
    "additional_errors": PDFset("LUXqed17_plus_PDF4LHC15_nnlo_100"),
    "luxseed": 123456789,
    "eps_base": 1e-2,
}

def test_parameters_init():

    fiatlux_runcard = FIATLUX_RUNCARD.copy()

    # we are not testing the photon here so we make it faster
    fiatlux_runcard['eps_base'] = 1e-1

    photon = Photon(TEST_THEORY, fiatlux_runcard, [1, 2, 3])
    alpha = Alpha(TEST_THEORY.get_description())

    np.testing.assert_equal(photon.replicas, [1, 2, 3])
    np.testing.assert_equal(photon.luxpdfset._name, FIATLUX_RUNCARD["luxset"].name)
    np.testing.assert_equal(photon.additional_errors._name, "LUXqed17_plus_PDF4LHC15_nnlo_100")
    np.testing.assert_equal(photon.luxseed, FIATLUX_RUNCARD["luxseed"])
    np.testing.assert_equal(photon.path_to_eko_photon, TEST_THEORY.path / "eko_photon.tar")
    np.testing.assert_equal(photon.q_in, 100.)
    np.testing.assert_almost_equal(
        alpha.alpha_em_ref, TEST_THEORY.get_description()["alphaqed"]
    )

def test_masses_init():
    alpha = Alpha(TEST_THEORY.get_description())
    np.testing.assert_equal(alpha.thresh_t, np.inf)
    np.testing.assert_almost_equal(alpha.thresh_b, 4.92)
    np.testing.assert_almost_equal(alpha.thresh_c, 1.51)

def test_set_thresholds_alpha_em():

    alpha = Alpha(TEST_THEORY.get_description())

    np.testing.assert_almost_equal(alpha.thresh[5], 91.2)
    np.testing.assert_almost_equal(alpha.thresh[4], 4.92)
    np.testing.assert_almost_equal(alpha.thresh[3], 51)
    np.testing.assert_almost_equal(alpha.alpha_thresh[5], 0.01)
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[4], alpha.alpha_em_fixed_flavor(4.92, 0.01, 91.2, 5)
    )
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[3],
        alpha.alpha_em_fixed_flavor(1.51, alpha.alpha_thresh[4], 4.92, 4),
    )
    np.testing.assert_equal(len(alpha.alpha_thresh), 3)
    np.testing.assert_equal(len(alpha.thresh), 3)

def test_betas():
    alpha = Alpha(TEST_THEORY.get_description())
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