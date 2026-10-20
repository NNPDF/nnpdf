"""Dependency-light checks for the default 'calibrated' (local-band) perturbation mode.

Run with ``python -m unittest discover -s validphys2/src/validphys/shapley/tests``.
Only NumPy is required; the band comes from a synthetic Stats class whose rows
are (central, lo, hi), as in test_consistency.py.
"""

import importlib.util
import os
from pathlib import Path
import sys
import types
import unittest

import numpy as np
from numpy.testing import assert_allclose


SOURCE = Path(os.environ.get('SHAPLEY_SOURCE_DIR', Path(__file__).resolve().parents[1]))
PACKAGE = '_shapley_local_band_tests'
package = types.ModuleType(PACKAGE)
package.__path__ = [str(SOURCE)]
sys.modules[PACKAGE] = package
spec = importlib.util.spec_from_file_location(f'{PACKAGE}.perturbation', SOURCE / 'perturbation.py')
perturbation = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = perturbation
spec.loader.exec_module(perturbation)


class SyntheticStats:
    def __init__(self, data):
        self.data = np.asarray(data)

    def central_value(self):
        return self.data[0]

    def errorbar68(self):
        return self.data[1], self.data[2]


MU, SIGMA = 0.3, 0.25
# Steeply falling PDF with a band that shrinks towards x = 1, and mu on a node.
X = np.unique(np.concatenate([np.logspace(-3, np.log10(0.95), 60), [MU]]))
CENTRAL = 2.0 * (1 - X) ** 5
HI = CENTRAL + 0.05 * CENTRAL * (1 + X)
LO = CENTRAL - 0.03 * CENTRAL * (1 + 2 * X)
CALIB = np.stack([CENTRAL, LO, HI])[:, None, :]          # (3 rows, 1 flavour, nx)
GV = np.repeat(CENTRAL[None, None, :], 2, axis=0)         # two members, one flavour


def perturb(mode, sign, cache=None, gv=GV, calib=CALIB):
    return perturbation.apply_gaussian_perturbation(
        gv, [0], MU, SIGMA, 1.0, X, mode=mode, xspace='logx',
        flavor_signs=np.full((gv.shape[0], 1), float(sign)),
        calibration_gv=calib, calibration_stats=SyntheticStats, amplitude_cache=cache,
    )[0, 0] - gv[0, 0]


class LocalBandTests(unittest.TestCase):
    def test_bump_stays_inside_local_band(self):
        up, down = perturb('calibrated', +1), perturb('calibrated', -1)
        self.assertTrue(np.all(up <= (HI - CENTRAL) * (1 + 1e-12)))
        self.assertTrue(np.all(down >= -(CENTRAL - LO) * (1 + 1e-12)))
        self.assertTrue(np.all(up >= 0) and np.all(down <= 0))

    def test_band_is_reached_only_at_the_centre(self):
        i = int(np.argmin(np.abs(X - MU)))
        up = perturb('calibrated', +1)
        assert_allclose(up[i], HI[i] - CENTRAL[i])
        self.assertTrue(np.all(np.delete(up, i) < np.delete(HI - CENTRAL, i)))

    def test_matches_centre_calibration_at_the_centre_only(self):
        i = int(np.argmin(np.abs(X - MU)))
        local, centre = perturb('calibrated', +1), perturb('calibrated_centre', +1)
        assert_allclose(local[i], centre[i])
        # Above the centre the band shrinks, so the centre-calibrated tail is larger.
        self.assertTrue(np.all(centre[i + 1:] > local[i + 1:]))

    def test_cache_is_reused_and_kept_apart_from_centre_amplitudes(self):
        cache = {}
        first = perturb('calibrated', +1, cache)
        self.assertIn(('local_band', 0), cache)
        stored = cache[('local_band', 0)]
        assert_allclose(perturb('calibrated', +1, cache), first)
        self.assertIs(cache[('local_band', 0)], stored)
        perturb('calibrated_centre', +1, cache)
        self.assertTrue(any(k != ('local_band', 0) for k in cache))
        assert_allclose(perturb('calibrated', +1, cache), first)

    def test_positivity_clip_keeps_the_bound(self):
        # Lower band edge below zero: the down shift would go negative and is clipped to 0.
        lo = CENTRAL - 1.5 * CENTRAL
        calib = np.stack([CENTRAL, lo, HI])[:, None, :]
        down = perturb('calibrated', -1, calib=calib)
        self.assertTrue(np.all(GV[0, 0] + down >= 0))
        self.assertTrue(np.all(down >= -(CENTRAL - lo) * (1 + 1e-12)))

    def multi(self, specs, signs):
        return perturbation.apply_multi_gaussian_perturbation(
            GV, specs, SIGMA, 1.0, X, mode='calibrated', xspace='logx',
            flavor_signs=np.tile(np.asarray(signs, dtype=float), (GV.shape[0], 1)),
            calibration_gv=CALIB, calibration_stats=SyntheticStats,
        )[0, 0] - GV[0, 0]

    def test_multi_bump_single_player_matches_single_x(self):
        for sign in (+1, -1):
            assert_allclose(self.multi([(0, MU, 0)], [sign]), perturb('calibrated', sign))

    def test_multi_bump_same_flavour_players_superpose(self):
        mu2 = 0.05
        both = self.multi([(0, MU, 0), (0, mu2, 1)], [+1, +1])
        g1 = perturbation.gaussian_profile(X, MU, SIGMA, 1.0, 'logx')
        g2 = perturbation.gaussian_profile(X, mu2, SIGMA, 1.0, 'logx')
        assert_allclose(both, (HI - CENTRAL) * (g1 + g2), atol=1e-14)


if __name__ == '__main__':
    unittest.main()
