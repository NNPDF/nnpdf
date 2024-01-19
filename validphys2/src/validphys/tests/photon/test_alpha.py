import numpy as np

from eko import beta
from eko.couplings import Couplings, expanded_qed
from eko.quantities.couplings import CouplingEvolutionMethod, CouplingsInfo
from eko.quantities.heavy_quarks import QuarkMassScheme
from validphys.api import API
from validphys.photon import constants
from validphys.photon.alpha import Alpha

from ..conftest import THEORY_QED


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

        # test all the thresholds
        np.testing.assert_almost_equal(alpha.alpha_em_ref, theory["alphaqed"])
        np.testing.assert_almost_equal(alpha.thresh[(5, 3)], theory["Qedref"])
        np.testing.assert_almost_equal(alpha.thresh[(4, 3)], theory["mb"])
        np.testing.assert_almost_equal(alpha.thresh[(4, 2)], constants.MTAU)
        np.testing.assert_almost_equal(alpha.thresh[(3, 2)], theory["mc"])
        np.testing.assert_almost_equal(alpha.thresh[(0, 2)], constants.MQL)
        np.testing.assert_almost_equal(alpha.thresh[(0, 1)], constants.MMU)
        np.testing.assert_almost_equal(alpha.thresh[(0, 0)], constants.ME)

        # test some values of alpha at the threshold points
        np.testing.assert_almost_equal(alpha.alphaem_thresh[(5, 3)], theory["alphaqed"])
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(4, 3)],
            alpha.alphaem_fixed_flavor(theory["mb"], alpha_ref, theory["Qedref"], 5, 3),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(4, 2)],
            alpha.alphaem_fixed_flavor(
                constants.MTAU, alpha.alphaem_thresh[(4, 3)], theory["mb"], 4, 3
            ),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(3, 2)],
            alpha.alphaem_fixed_flavor(
                theory["mc"], alpha.alphaem_thresh[(4, 2)], constants.MTAU, 4, 2
            ),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(0, 2)],
            alpha.alphaem_fixed_flavor(
                constants.MQL, alpha.alphaem_thresh[(3, 2)], theory["mc"], 3, 2
            ),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(0, 1)],
            alpha.alphaem_fixed_flavor(
                constants.MMU, alpha.alphaem_thresh[(0, 2)], constants.MQL, 0, 2
            ),
        )
        np.testing.assert_almost_equal(
            alpha.alphaem_thresh[(0, 0)],
            alpha.alphaem_fixed_flavor(
                constants.ME, alpha.alphaem_thresh[(0, 1)], constants.MMU, 0, 1
            ),
        )

        np.testing.assert_equal(len(alpha.alphaem_thresh), 7)
        np.testing.assert_equal(len(alpha.thresh), 7)


def test_couplings_exa():
    """
    Benchmark the exact running of alpha:
    Testing both the FFS and the VFS with the results from EKO.
    Test also the values of the couplings at the threshold values.
    """
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()
    for k in [1]:
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
                alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5, 3),
                alpha.alpha_em(q),
                rtol=5e-6,
            )
            np.testing.assert_allclose(
                alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5, 3),
                eko_alpha.compute_exact_alphaem_running(
                    np.array([0.118, alpha_ref]) / (4 * np.pi), 5, 3, theory["Qref"] ** 2, q**2
                )[1]
                * 4
                * np.pi,
                rtol=1e-7,
            )

        # TODO: these tests have to be switched on once the varying nl is implemented in eko

        for q in [1, 2, 3, 4]:
            np.testing.assert_allclose(
                alpha.alpha_em(q), eko_alpha.a_em(q**2) * 4 * np.pi, rtol=5e-6
            )
        for nf in range(3, theory["MaxNfAs"]):
            np.testing.assert_allclose(
                alpha.alphaem_thresh[(nf, 2 if nf == 3 else 3)],
                eko_alpha.a_em(mass_list[nf - 3] ** 2, nf) * 4 * np.pi,
                rtol=3e-7,
            )


def test_exa_interpolation():
    """test the accuracy of the alphaem interpolation"""
    test_theory = API.theoryid(theoryid=THEORY_QED)
    theory = test_theory.get_description()

    alpha = Alpha(theory, 1e8)
    for q in np.geomspace(5.0, 1e4, 1000, endpoint=True):
        np.testing.assert_allclose(
            alpha.alpha_em(q),
            alpha.alphaem_fixed_flavor(q, theory["alphaqed"], theory["Qref"], 5, 3),
            rtol=1e-5,
        )
    for q in np.geomspace(1.778, 4.9, 20):
        np.testing.assert_allclose(
            alpha.alpha_em(q),
            alpha.alphaem_fixed_flavor(q, alpha.alphaem_thresh[(4, 3)], alpha.thresh_b, 4, 3),
            rtol=1e-5,
        )
    for q in np.geomspace(0.2, 0.8, 20):
        np.testing.assert_allclose(
            alpha.alpha_em(q),
            alpha.alphaem_fixed_flavor(q, alpha.alphaem_thresh[(0, 2)], constants.MQL, 0, 2),
            rtol=1e-5,
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
            alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5, 3),
            alpha.alpha_em(q),
            rtol=1e-10,
        )
        np.testing.assert_allclose(
            alpha.alphaem_fixed_flavor(q, alpha_ref, theory["Qref"], 5, 3),
            expanded_qed(
                alpha_ref / (4 * np.pi),
                theory["QED"],
                beta.beta_qed((0, 2), 5, 3),
                [beta.b_qed((0, i + 2), 5, 3) for i in range(theory["QED"])],
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
        -0.0,
        -0.10610329539459688,
        -0.21220659078919377,
        -0.42441318157838753,
        -0.5658842421045168,
        -0.6719875374991137,
        -0.7073553026306458,
    ]
    vec_beta1 = [
        -0.0,
        -0.025330295910584444,
        -0.05066059182116889,
        -0.06754745576155852,
        -0.0825580014863493,
        -0.10788829739693376,
        -0.10882645650473316,
    ]
    for i, (nf, nl) in enumerate([(0, 0), (0, 1), (0, 2), (3, 2), (4, 2), (4, 3), (5, 3)]):
        np.testing.assert_allclose(alpha.betas_qed[(nf, nl)][0], vec_beta0[i], rtol=1e-7)
        np.testing.assert_allclose(alpha.betas_qed[(nf, nl)][1], vec_beta1[i], rtol=1e-7)
