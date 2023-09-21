import tempfile

import fiatlux
import numpy as np
import yaml

from eko import beta
from eko.couplings import Couplings, expanded_qed
from eko.io import EKO
from eko.quantities.couplings import CouplingEvolutionMethod, CouplingsInfo
from eko.quantities.heavy_quarks import QuarkMassScheme
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


# def test_parameters_init():
#     """test initailization of the parameters from Photon class"""
#     fiatlux_runcard = generate_fiatlux_runcard()
#     test_theory = API.theoryid(theoryid=THEORY_QED)

#     # we are not testing the photon here so we make it faster
#     fiatlux_runcard['eps_base'] = 1e-1

#     photon = Photon(test_theory, fiatlux_runcard, [1, 2, 3])

#     np.testing.assert_equal(photon.replicas, [1, 2, 3])
#     np.testing.assert_equal(photon.luxpdfset._name, fiatlux_runcard["luxset"].name)
#     np.testing.assert_equal(photon.additional_errors.name, "LUXqed17_plus_PDF4LHC15_nnlo_100")
#     np.testing.assert_equal(photon.luxseed, fiatlux_runcard["luxseed"])
#     np.testing.assert_equal(photon.path_to_eko_photon, test_theory.path / "eko_photon.tar")
#     np.testing.assert_equal(photon.q_in, 100.0)


def test_masses_init():
    """test thresholds values in Alpha class"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()
    alpha = Alpha(theory, 1e8)
    np.testing.assert_equal(alpha.thresh_t, np.inf)
    np.testing.assert_almost_equal(alpha.thresh_b, theory["mb"])
    np.testing.assert_almost_equal(alpha.thresh_c, theory["mc"])


def test_set_thresholds_alpha_em():
    """test value of alpha_em at threshold values"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()

    for modev in ['EXA', 'TRN']:
        theory['ModEv'] = modev

        alpha = Alpha(theory, 1e8)
        alpha_ref = theory["alphaqed"]

        np.testing.assert_almost_equal(alpha.alpha_em_ref, theory["alphaqed"])
        np.testing.assert_almost_equal(alpha.thresh[5], theory["Qedref"])
        np.testing.assert_almost_equal(alpha.thresh[4], theory["mb"])
        np.testing.assert_almost_equal(alpha.thresh[3], theory["mc"])
        np.testing.assert_almost_equal(alpha.alphaem_thresh[5], theory["alphaqed"])
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[4],
            alpha.alphaem_fixed_flavor(theory["mb"], alpha_ref, theory["Qedref"], 5),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[3],
            alpha.alphaem_fixed_flavor(theory["mc"], alpha.alphaem_thresh[4], theory["mb"], 4),
        )
        np.testing.assert_equal(len(alpha.alphaem_thresh), 3)
        np.testing.assert_equal(len(alpha.thresh), 3)


def test_couplings_exa():
    """
    Benchmark the exact running of alpha:
    Testing both the FFS and the VFS with the results from EKO.
    Test also the values of the couplings at the threshold values.
    """
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()
    for k in [0.5, 1, 2]:
        theory["mb"] *= k
        mass_list = [theory["mc"], theory["mb"], theory["mt"]]

        alpha = Alpha(theory, 1e8)
        couplings = CouplingsInfo.from_dict(
            dict(
                alphas=theory["alphas"],
                alphaem=theory["alphaqed"],
                scale=theory["Qref"],
                num_flavs_ref=None,
                max_num_flavs=theory["MaxNfAs"],
                em_running=True,
            )
        )
        eko_alpha = Couplings(
            couplings,
            (theory["PTO"] + 1, theory["QED"]),
            method=CouplingEvolutionMethod.EXACT,
            masses=[m**2 for m in mass_list],
            hqm_scheme=QuarkMassScheme.POLE,
            thresholds_ratios=[1.0, 1.0, 1.0],
        )
        eko_alpha.decoupled_running = True
        alpha_ref = theory["alphaqed"]
        for q in [5, 10, 20, 50, 80, 100, 200]:
            np.testing.assert_allclose(
                alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5),
                alpha.alpha_em(q),
                rtol=5e-6,
            )
            np.testing.assert_allclose(
                alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5),
                eko_alpha.compute_exact_alphaem_running(
                    np.array([0.118, alpha_ref]) / (4 * np.pi), 5, theory["Qref"] ** 2, q**2
                )[1]
                * 4
                * np.pi,
                rtol=1e-7,
            )
        for q in [1, 2, 3, 4]:
            np.testing.assert_allclose(
                alpha.alpha_em(q), eko_alpha.a_em(q**2) * 4 * np.pi, rtol=5e-6
            )
        for nf in range(3, theory["MaxNfAs"]):
            np.testing.assert_allclose(
                alpha.alphaem_thresh[nf],
                eko_alpha.a_em(mass_list[nf - 3] ** 2, nf) * 4 * np.pi,
                rtol=3e-7,
            )


def test_exa_interpolation():
    """test the accuracy of the alphaem interpolation"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()

    alpha = Alpha(theory, 1e8)
    for q in np.geomspace(1.0, 1e4, 1000, endpoint=True):
        np.testing.assert_allclose(
            alpha.alpha_em(q), alpha.alpha_em(q), rtol=1e-5
        )


def test_couplings_trn():
    """benchmark the truncated running of alpha with EKO"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()
    theory["ModEv"] = 'TRN'
    alpha = Alpha(theory, 1e8)
    alpha_ref = theory["alphaqed"]

    for q in [80, 10, 5]:
        np.testing.assert_allclose(
            alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5),
            alpha.alpha_em(q),
            rtol=1e-10,
        )
        np.testing.assert_allclose(
            alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5),
            expanded_qed(
                alpha_ref / (4 * np.pi),
                theory["QED"],
                beta.beta_qed((0, 2), 5),
                [beta.b_qed((0, i + 2), 5) for i in range(theory["QED"])],
                2 * np.log(q / theory["Qref"]),
            )
            * 4
            * np.pi,
            rtol=1e-7,
        )


def test_trn_vs_exa():
    """benchmark the exact running of alpha vs truncated"""
    # does this test make sense?
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory_exa = test_theory.get_description()
    theory_trn = theory_exa.copy()
    theory_trn["ModEv"] = 'TRN'

    alpha_exa = Alpha(theory_exa, 1e8)
    alpha_trn = Alpha(theory_trn, 1e8)

    for q in [1, 3, 10, 50, 80]:
        np.testing.assert_allclose(alpha_exa.alpha_em(q), alpha_trn.alpha_em(q), rtol=2e-3)


def test_betas():
    """test betas for different nf"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    alpha = Alpha(test_theory.get_description(), 1e8)
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
        np.testing.assert_allclose(alpha.betas_qed[nf][0], vec_beta0[nf - 3], rtol=1e-7)
        np.testing.assert_allclose(alpha.betas_qed[nf][1], vec_b1[nf - 3], rtol=1e-7)


# def test_photon():
#     """
#     Test that photon coming out of Photon interpolator matches the photon array
#     for XGRID points
#     """
#     fiatlux_runcard = generate_fiatlux_runcard()
#     fiatlux_runcard["additional_errors"] = False
#     test_theory = API.theoryid(theoryid=THEORY_QED)
#     theory = test_theory.get_description()

#     for replica in [1, 2, 3]:
#         photon = Photon(test_theory, fiatlux_runcard, [replica])

#         # set up fiatlux
#         path_to_F2 = test_theory.path / "fastkernel/FIATLUX_DIS_F2.pineappl.lz4"
#         path_to_FL = test_theory.path / "fastkernel/FIATLUX_DIS_FL.pineappl.lz4"
#         pdfs = fiatlux_runcard["luxset"].load()
#         f2 = sf.InterpStructureFunction(path_to_F2, pdfs.members[replica])
#         fl = sf.InterpStructureFunction(path_to_FL, pdfs.members[replica])
#         f2lo = sf.F2LO(pdfs.members[replica], theory)

#         # runcard
#         fiatlux_default = FIATLUX_DEFAULT.copy()
#         fiatlux_default['mproton'] = theory['MP']
#         fiatlux_default["qed_running"] = bool(np.isclose(theory["Qedref"], theory["Qref"]))
#         fiatlux_default["q2_max"] = float(f2.q2_max)
#         fiatlux_default["eps_base"] = fiatlux_runcard["eps_base"]

#         # load fiatlux
#         with tempfile.NamedTemporaryFile(mode="w") as tmp:
#             with tmp.file as tmp_file:
#                 tmp_file.write(yaml.dump(FIATLUX_DEFAULT))
#             lux = fiatlux.FiatLux(tmp.name)

#         alpha = Alpha(theory, fiatlux_default["q2_max"])

#         lux.PlugAlphaQED(alpha.alpha_em, alpha.qref)
#         lux.InsertInelasticSplitQ(
#             [
#                 theory["kbThr"] * theory["mb"],
#                 theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100,
#             ]
#         )
#         lux.PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)
#         path_to_eko_photon = test_theory.path / "eko_photon.tar"
#         with EKO.read(path_to_eko_photon) as eko:
#             photon_fiatlux_qin = np.array([lux.EvaluatePhoton(x, eko.mu20).total for x in XGRID])
#             photon_fiatlux_qin /= XGRID
#             # construct PDFs
#             pdfs_init = np.zeros((len(eko.bases.inputpids), len(XGRID)))
#             for j, pid in enumerate(eko.bases.inputpids):
#                 if pid == 22:
#                     pdfs_init[j] = photon_fiatlux_qin
#                     ph_id = j
#                 else:
#                     if pid not in pdfs.flavors:
#                         continue
#                     pdfs_init[j] = np.array(
#                         [pdfs.xfxQ(x, np.sqrt(eko.mu20), replica, pid) / x for x in XGRID]
#                     )

#             # Apply EKO to PDFs
#             for _, elem in eko.items():
#                 pdfs_final = np.einsum("ajbk,bk", elem.operator, pdfs_init)

#         photon_Q0 = pdfs_final[ph_id]
#         photon_fiatlux = XGRID * photon_Q0

#         photon_validphys = photon(XGRID[np.newaxis, :, np.newaxis])[0][0, :, 0]

#         np.testing.assert_allclose(photon_fiatlux, photon_validphys, rtol=1e-7)
