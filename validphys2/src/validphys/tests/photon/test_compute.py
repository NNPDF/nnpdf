import tempfile

import fiatlux
import numpy as np
import yaml

from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.api import API
from validphys.core import PDF as PDFset
from validphys.loader import FallbackLoader
from validphys.photon import structure_functions as sf
from validphys.photon.compute import FIATLUX_DEFAULT, Alpha, Photon

from ..conftest import PDF, THEORY_QED


def generate_fiatlux_runcard():
    return {
        "luxset": PDFset(PDF),
        # check if "LUXqed17_plus_PDF4LHC15_nnlo_100" is installed
        "additional_errors": FallbackLoader().check_pdf("LUXqed17_plus_PDF4LHC15_nnlo_100"),
        "luxseed": 123456789,
        "eps_base": 1e-2,  # using low precision to speed up tests
    }


def test_parameters_init():
    "test initailization of the parameters from Photon class"
    fiatlux_runcard = generate_fiatlux_runcard()
    test_theory = API.theoryid(theoryid=THEORY_QED)

    # we are not testing the photon here so we make it faster
    fiatlux_runcard['eps_base'] = 1e-1

    photon = Photon(test_theory, fiatlux_runcard, [1, 2, 3])

    np.testing.assert_equal(photon.replicas, [1, 2, 3])
    np.testing.assert_equal(photon.luxpdfset._name, fiatlux_runcard["luxset"].name)
    np.testing.assert_equal(photon.additional_errors.name, "LUXqed17_plus_PDF4LHC15_nnlo_100")
    np.testing.assert_equal(photon.luxseed, fiatlux_runcard["luxseed"])
    np.testing.assert_equal(photon.path_to_eko_photon, test_theory.path / "eko_photon.tar")
    np.testing.assert_equal(photon.q_in, 100.0)


def test_masses_init():
    "test thresholds values in Alpha class"
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()
    alpha = Alpha(theory)
    np.testing.assert_equal(alpha.thresh_t, np.inf)
    np.testing.assert_almost_equal(alpha.thresh_b, theory["mb"])
    np.testing.assert_almost_equal(alpha.thresh_c, theory["mc"])


def test_set_thresholds_alpha_em():
    "test value of alpha_em at threshold values"
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()

    alpha = Alpha(theory)

    np.testing.assert_almost_equal(alpha.alpha_em_ref, theory["alphaqed"])
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
    "test betas for different nf"
    test_theory = API.theoryid(theoryid=THEORY_QED)
    alpha = Alpha(test_theory.get_description())
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
    """test that photon coming out of Photon interpolator matches the photon array
    for XGRID points
    """
    fiatlux_runcard = generate_fiatlux_runcard()
    fiatlux_runcard["additional_errors"] = False
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()

    for replica in [1, 2, 3]:
        photon = Photon(test_theory, fiatlux_runcard, [replica])

        # set up fiatlux
        path_to_F2 = test_theory.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
        path_to_FL = test_theory.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"
        pdfs = fiatlux_runcard["luxset"].load()
        f2 = sf.InterpStructureFunction(path_to_F2, pdfs.members[replica])
        fl = sf.InterpStructureFunction(path_to_FL, pdfs.members[replica])
        f2lo = sf.F2LO(pdfs.members[replica], theory)

        # runcard
        fiatlux_default = FIATLUX_DEFAULT.copy()
        fiatlux_default['mproton'] = theory['MP']
        fiatlux_default["qed_running"] = bool(np.isclose(theory["Qedref"], theory["Qref"]))
        fiatlux_default["q2_max"] = float(f2.q2_max)
        fiatlux_default["eps_base"] = fiatlux_runcard["eps_base"]

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
        path_to_eko_photon = test_theory.path / "eko_photon.tar"
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
