"""Centre calibration for the dev-only paths: amplitude cache, symmetric mode, multi-player.

Reuses the synthetic set-up of test_consistency. Run with
``python -m unittest discover -s validphys2/src/validphys/shapley/tests``.
"""

import types
import unittest

import numpy as np
from numpy.testing import assert_allclose

from test_consistency import SyntheticObservable, SyntheticStats, analyzer_module


class SyntheticStatsStd(SyntheticStats):
    def std_error(self):
        return 0.5 * (self.data[2] - self.data[1])


FLAVOR_INFO = dict(names=['valence', 'test'], indices=[[3, 6], 9], n_flavors=2)
MODES = ('calibrated_centre', 'calibrated_symmetric')


def observables(hadronic=False):
    return [SyntheticObservable('a', [.1, .3], hadronic),
            SyntheticObservable('b', [.2, .3], hadronic)]


def pdf():
    return types.SimpleNamespace(stats_class=SyntheticStatsStd, error_type='replicas')


def scalar(basis, enforce, cache, hadronic=False):
    return analyzer_module.NNPDFShapleyAnalyzer(
        pdf(), observables(hadronic), FLAVOR_INFO, basis=basis,
        enforce_sumrules=enforce, member_mode='central', cache_calibration=cache,
    )


def multiplayer(basis, enforce, cache, mode, players, hadronic=False):
    return analyzer_module.NNPDFShapleyAnalyzerMultiPlayer(
        pdf(), observables(hadronic), FLAVOR_INFO, players, .1, .2, mode, 'linear',
        basis=basis, enforce_sumrules=enforce, member_mode='central',
        cache_calibration=cache,
    )


SIGNS = (np.array([[1., -1.]]), np.array([[-1., 1.]]))


class CacheTests(unittest.TestCase):
    def test_cache_is_bit_identical_scalar(self):
        for basis in ('flavor', 'evolution'):
            for enforce in (False, True):
                for mode in MODES:
                    off, on = scalar(basis, enforce, False), scalar(basis, enforce, True)
                    for coalition in ([0], [1], [0, 1]):
                        for signs in SIGNS:
                            with self.subTest(basis=basis, enforce=enforce, mode=mode,
                                              coalition=coalition):
                                kwargs = dict(mode=mode, random_sign_matrix=signs)
                                # Evaluate twice with the cache on so the second call hits it.
                                on._evaluate_chi2(coalition, .15, .1, .2, **kwargs)
                                assert_allclose(
                                    on._evaluate_chi2(coalition, .15, .1, .2, **kwargs),
                                    off._evaluate_chi2(coalition, .15, .1, .2, **kwargs),
                                    rtol=0, atol=0,
                                )

    def test_cache_is_bit_identical_multiplayer(self):
        players = [(0, .15), (1, .25)]
        for basis in ('flavor', 'evolution'):
            for mode in MODES:
                off = multiplayer(basis, True, False, mode, players)
                on = multiplayer(basis, True, True, mode, players)
                for coalition in ([0], [1], [0, 1]):
                    for signs in SIGNS:
                        with self.subTest(basis=basis, mode=mode, coalition=coalition):
                            on._evaluate_chi2(coalition, None, None, None,
                                              random_sign_matrix=signs)
                            assert_allclose(
                                on._evaluate_chi2(coalition, None, None, None,
                                                  random_sign_matrix=signs),
                                off._evaluate_chi2(coalition, None, None, None,
                                                   random_sign_matrix=signs),
                                rtol=0, atol=0,
                            )


class GridIndependenceTests(unittest.TestCase):
    def test_symmetric_mode_same_pdf_on_different_fk_grids(self):
        for basis in ('flavor', 'evolution'):
            for enforce in (False, True):
                for hadronic in (False, True):
                    with self.subTest(basis=basis, enforce=enforce, hadronic=hadronic):
                        analyzer = scalar(basis, enforce, True, hadronic)
                        analyzer._evaluate_chi2([0, 1], .15, .1, .2,
                                                mode='calibrated_symmetric')
                        a, b = analyzer.observables
                        assert_allclose(a.last_grid[:, :, -1], b.last_grid[:, :, -1],
                                        rtol=0, atol=1e-14)

    def test_multiplayer_same_pdf_on_different_fk_grids(self):
        # Players at two different centres, neither on a grid node.
        players = [(0, .15), (1, .25)]
        for basis in ('flavor', 'evolution'):
            for enforce in (False, True):
                for hadronic in (False, True):
                    for mode in MODES:
                        with self.subTest(basis=basis, enforce=enforce,
                                          hadronic=hadronic, mode=mode):
                            analyzer = multiplayer(basis, enforce, True, mode, players,
                                                   hadronic)
                            analyzer._evaluate_chi2([0, 1], None, None, None,
                                                    random_sign_matrix=SIGNS[0])
                            a, b = analyzer.observables
                            assert_allclose(a.last_grid[:, :, -1], b.last_grid[:, :, -1],
                                            rtol=0, atol=1e-14)


if __name__ == '__main__':
    unittest.main()
