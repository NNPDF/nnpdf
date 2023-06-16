import tempfile

import fiatlux
import numpy as np
import yaml

from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.api import API
from validphys.core import PDF as PDFset
from validphys.photon import structure_functions as sf
from validphys.photon.compute import Alpha, Photon

from ..conftest import PDF

TEST_THEORY = API.theoryid(theoryid=398)

FIATLUX_RUNCARD = {
    "luxset": PDFset(PDF),
    "additional_errors": PDFset("LUXqed17_plus_PDF4LHC15_nnlo_100"),
    "luxseed": 123456789,
    "eps_base": 1e-2, # using low precision to speed up tests
}

FIATLUX_DEFAULT = {
    "apfel": False,
    "eps_rel": 1e-1,  # extra precision on any single integration.
    "mum_proton": 2.792847356,  # proton magnetic moment, from
    # http://pdglive.lbl.gov/DataBlock.action?node=S016MM which itself
    # gets it from arXiv:1203.5425 (CODATA)
    # the elastic param type, options:
    # dipole
    # A1_world_spline
    # A1_world_pol_spline
    "elastic_param": "A1_world_pol_spline",
    "elastic_electric_rescale": 1,
    "elastic_magnetic_rescale": 1,
    # the inelastic param type, options:
    "inelastic_param": "LHAPDF_Hermes_ALLM_CLAS",  # Hermes_ALLM_CLAS, LHAPDF_Hermes_ALLM_CLAS
    "rescale_r_twist4": 0,
    "rescale_r": 1,
    "allm_limits": 0,
    "rescale_non_resonance": 1,
    "rescale_resonance": 1,
    "use_mu2_as_upper_limit": False,
    "q2min_inel_override": 0.0,
    "q2max_inel_override": 1e300,
    "lhapdf_transition_q2": 9,
    # general
    "verbose": False,
}


def test_parameters_init():
    fiatlux_runcard = FIATLUX_RUNCARD.copy()

    # we are not testing the photon here so we make it faster
    fiatlux_runcard['eps_base'] = 1e-1

    photon = Photon(TEST_THEORY, fiatlux_runcard, [1, 2, 3])
    alpha = Alpha(TEST_THEORY.get_description())

    np.testing.assert_equal(photon.replicas, [1, 2, 3])
    np.testing.assert_equal(photon.luxpdfset._name, FIATLUX_RUNCARD["luxset"].name)
    np.testing.assert_equal(photon.additional_errors.name, "LUXqed17_plus_PDF4LHC15_nnlo_100")
    np.testing.assert_equal(photon.luxseed, FIATLUX_RUNCARD["luxseed"])
    np.testing.assert_equal(photon.path_to_eko_photon, TEST_THEORY.path / "eko_photon.tar")
    np.testing.assert_equal(photon.q_in, 100.0)
    np.testing.assert_almost_equal(alpha.alpha_em_ref, TEST_THEORY.get_description()["alphaqed"])


def test_masses_init():
    theory = TEST_THEORY.get_description()
    alpha = Alpha(theory)
    np.testing.assert_equal(alpha.thresh_t, np.inf)
    np.testing.assert_almost_equal(alpha.thresh_b, theory["mb"])
    np.testing.assert_almost_equal(alpha.thresh_c, theory["mc"])


def test_set_thresholds_alpha_em():
    theory = TEST_THEORY.get_description()

    alpha = Alpha(theory)

    np.testing.assert_almost_equal(alpha.thresh[5], theory["Qedref"])
    np.testing.assert_almost_equal(alpha.thresh[4], theory["mb"])
    np.testing.assert_almost_equal(alpha.thresh[3], theory["mc"])
    np.testing.assert_almost_equal(alpha.alpha_thresh[5], theory["alphaqed"])
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[4],
        alpha.alpha_em_fixed_flavor(theory["mb"], theory["alphaqed"], theory["Qedref"], 5),
    )
    np.testing.assert_almost_equal(
        alpha.alpha_thresh[3],
        alpha.alpha_em_fixed_flavor(theory["mc"], alpha.alpha_thresh[4], theory["mb"], 4),
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


def test_photon():
    fiatlux_runcard = FIATLUX_RUNCARD.copy()
    fiatlux_runcard["additional_errors"] = False
    theory = TEST_THEORY.get_description()

    for replica in [1, 2, 3]:
        photon = Photon(TEST_THEORY, fiatlux_runcard, [replica])

        # set up fiatlux
        path_to_F2 = TEST_THEORY.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
        path_to_FL = TEST_THEORY.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"
        pdfs = FIATLUX_RUNCARD["luxset"].load()
        f2 = sf.InterpStructureFunction(path_to_F2, pdfs.members[replica])
        fl = sf.InterpStructureFunction(path_to_FL, pdfs.members[replica])
        f2lo = sf.F2LO(pdfs.members[replica], theory)

        # runcard
        fiatlux_default = FIATLUX_DEFAULT.copy()
        fiatlux_default['mproton'] = theory['MP']
        fiatlux_default["qed_running"] = bool(np.isclose(theory["Qedref"], theory["Qref"]))
        fiatlux_default["q2_max"] = float(f2.q2_max)
        fiatlux_default["eps_base"] = FIATLUX_RUNCARD["eps_base"]

        # load fiatlux
        with tempfile.NamedTemporaryFile(mode="w") as tmp:
            with tmp.file as tmp_file:
                tmp_file.write(yaml.dump(FIATLUX_DEFAULT))
            lux = fiatlux.FiatLux(tmp.name)

        alpha = Alpha(theory)

        lux.PlugAlphaQED(alpha.alpha_em, alpha.qref)
        lux.InsertInelasticSplitQ(
            [
                theory["kbThr"] * theory["mb"],
                theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100,
            ]
        )
        lux.PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)
        path_to_eko_photon = TEST_THEORY.path / "eko_photon.tar"
        with EKO.read(path_to_eko_photon) as eko:
            photon_fiatlux_qin = np.array([lux.EvaluatePhoton(x, eko.mu20).total for x in XGRID])
            photon_fiatlux_qin /= XGRID
            # construct PDFs
            pdfs_init = np.zeros((len(eko.bases.inputpids), len(XGRID)))
            for j, pid in enumerate(eko.bases.inputpids):
                if pid == 22:
                    pdfs_init[j] = photon_fiatlux_qin
                    ph_id = j
                else:
                    if pid not in pdfs.flavors:
                        continue
                    pdfs_init[j] = np.array(
                        [pdfs.xfxQ(x, np.sqrt(eko.mu20), replica, pid) / x for x in XGRID]
                    )

            # Apply EKO to PDFs
            for _, elem in eko.items():
                pdfs_final = np.einsum("ajbk,bk", elem.operator, pdfs_init)

        photon_Q0 = pdfs_final[ph_id]
        photon_fiatlux = XGRID * photon_Q0

        photon_validphys = photon(XGRID[np.newaxis, :, np.newaxis])[0][0, :, 0]

        np.testing.assert_allclose(photon_fiatlux, photon_validphys, rtol=1e-7)
