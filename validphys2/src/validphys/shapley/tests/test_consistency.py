"""Dependency-light regressions for the Shapley perturbation prescription.

Run with ``python -m unittest discover -s validphys2/src/validphys/shapley/tests``.
Only NumPy is required. PDF loading, rotations and convolution are synthetic;
the production perturbation, normalization, analyzer and coalition-sum code
are loaded unchanged. These tests do not replace an LHAPDF/FK integration run.
Set SHAPLEY_SOURCE_DIR to check another revision of the three source modules.
"""

import importlib.util
from itertools import combinations
import os
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose


SOURCE = Path(os.environ.get('SHAPLEY_SOURCE_DIR', Path(__file__).resolve().parents[1]))
PACKAGE = '_shapley_consistency_tests'
package = types.ModuleType(PACKAGE)
package.__path__ = [str(SOURCE)]
sys.modules[PACKAGE] = package


def load_module(name):
    spec = importlib.util.spec_from_file_location(f'{PACKAGE}.{name}', SOURCE / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


perturbation = load_module('perturbation')
sumrules = load_module('sumrules')


class SyntheticStats:
    def __init__(self, data):
        self.data = np.asarray(data)

    def central_value(self):
        return self.data[0]

    def errorbar68(self):
        return self.data[1], self.data[2]


def pdf_values(pdf, target, n_replicas=None, member_mode='all'):
    x = np.asarray(target.xgrid)
    central = np.full((14, len(x)), 5.0)
    central[0] = 0.01
    central[1] = 0.49
    central[2] = 0.5
    central[3:9] = np.array([3, 1, 3, 3, 3, 3])[:, None] * x
    # Vary the uncertainty with both x and flavour to expose wrong grid or
    # flavour indexing. Deliberately asymmetric confidence bounds.
    spread = 0.01 * (1 + np.arange(14))[:, None] * x
    members = np.stack((central, central - spread, central + 2 * spread))
    return members[:1] if member_mode == 'central' else members


def subset_pdf_values(pdf, target, **kwargs):
    return pdf_values(pdf, target, **kwargs)[:, target.flavor_indices]


def stub_module(name, **attributes):
    module = types.ModuleType(name)
    module.__dict__.update(attributes)
    return module


# Stub only external dependencies. All numerical code under review is real.
stubs = {
    'validphys.convolution': stub_module('validphys.convolution', FK_FLAVOURS=np.arange(14)),
    'validphys.pdfbases': stub_module(
        'validphys.pdfbases', ALL_FLAVOURS=np.arange(14),
        evolution=types.SimpleNamespace(
            _to_indexes=lambda _: np.arange(14), from_flavour_mat=np.eye(14),
        ),
    ),
    'shapley_values': stub_module('shapley_values', ExactShapley=None, plot_shapley_bar=None),
    'matplotlib': stub_module('matplotlib'),
    'matplotlib.pyplot': stub_module('matplotlib.pyplot'),
    f'{PACKAGE}.setup': stub_module(
        f'{PACKAGE}.setup', get_pdf_grid_values=subset_pdf_values,
        get_pdf_flavor_grid_values=pdf_values, get_pdf_grid_values_all14=pdf_values,
        FLAVOR_PDG_NAMES={},
    ),
}
with patch.dict(sys.modules, stubs):
    analyzer_module = load_module('analyzer')


class SyntheticObservable:
    def __init__(self, name, grid, hadronic=False):
        self.name, self.Q0, self.ndata = name, 1.65, 1
        self.fk_entries = [types.SimpleNamespace(
            Q0=self.Q0, xgrid=np.asarray(grid), hadronic=hadronic,
            flavor_indices=np.array([2, 3, 4, 6, 9]),
        )]

    def rotate_to_evolution(self, grids):
        # Identity synthetic rotation; selection still exercises local ordering.
        entry = self.fk_entries[0]
        return grids if entry.hadronic else [grids[0][:, entry.flavor_indices]]

    def chi2(self, grids):
        self.last_grid = grids[0].copy()
        return np.sum(grids[0][:, :, -1], axis=1) ** 2

    def chi2_from_flavor(self, grids):
        return self.chi2(self.rotate_to_evolution(grids))


def make_analyzer(basis='evolution', enforce=True, vector=False, hadronic=False):
    cls = (analyzer_module.NNPDFShapleyAnalyzerVecX if vector
           else analyzer_module.NNPDFShapleyAnalyzer)
    pdf = types.SimpleNamespace(stats_class=SyntheticStats, error_type='replicas')
    observables = [SyntheticObservable('a', [.1, .3], hadronic),
                   SyntheticObservable('b', [.2, .3], hadronic)]
    kwargs = dict(basis=basis, enforce_sumrules=enforce, member_mode='central')
    if vector:
        kwargs.update(x_values=[.15], vec_sigma=.1, vec_amplitude=.2,
                      vec_mode='calibrated', vec_xspace='linear')
    return cls(pdf, observables,
               dict(names=['valence', 'test'], indices=[[3, 6], 9], n_flavors=2),
               **kwargs)


class SumruleTests(unittest.TestCase):
    def test_signed_and_independent_valence_integrals(self):
        x, weights = sumrules.gen_integration_input()
        gv = np.zeros((2, 14, len(x)))
        gv[:, 0] = .01
        gv[:, 1] = .49
        gv[:, 2] = np.array([.5, -.5])[:, None]
        gv[:, 3:9] = np.array([[2, 4, 5, 6, 7, 8], [-2, -4, -5, -6, -7, -8]])[:, :, None] * x
        norm = sumrules.compute_sumrule_normalization(gv, x, weights)
        corrected = gv * norm[:, :, None]
        momentum = np.einsum('rfx,x->rf', corrected, weights)
        number = np.einsum('rfx,x,x->rf', corrected, 1 / x, weights)
        assert_allclose(momentum[:, :3].sum(axis=1), 1, atol=1e-12)
        assert_allclose(number[:, 3:9], np.tile([3, 1, 3, 3, 3, 3], (2, 1)), atol=1e-12)
        assert_allclose(sumrules.compute_sumrule_normalization(corrected, x, weights), 1,
                        atol=1e-12)

    def test_singular_normalization_is_explicit(self):
        x, weights = sumrules.gen_integration_input()
        with self.assertRaises(ValueError):
            sumrules.compute_sumrule_normalization(np.zeros((1, 14, len(x))), x, weights)


class PerturbationTests(unittest.TestCase):
    def test_explicit_centre_calibration_is_grid_independent(self):
        centre = np.array([[1.5], [1.2], [1.8]])
        for sign in (1., -1.):
            outputs = []
            for x in (np.array([.1, .3]), np.array([.2, .3])):
                outputs.append(perturbation.apply_gaussian_perturbation(
                    np.full((1, 1, 2), 10.), [0], .15, .1, 1, x,
                    mode='calibrated', flavor_signs=np.array([[sign]]),
                    calibration_at_mu=centre, calibration_stats=SyntheticStats,
                )[0, 0, -1])
            assert_allclose(outputs[0], outputs[1], rtol=0, atol=1e-14)

    def test_clipping_can_leave_a_linear_symmetrized_response(self):
        # X(f)=(f-2)^2, f0=0. Clipping gives f+=eps, f-=0, hence
        # [X(f+)+X(f-)]/2-X(0) = -2 eps + eps^2/2.
        eps = 1e-5
        values = []
        for sign in (1., -1.):
            value = perturbation.apply_gaussian_perturbation(
                np.zeros((1, 1, 1)), [0], .2, .1, eps, np.array([.2]),
                flavor_signs=np.array([[sign]]),
            )[0, 0, 0]
            values.append((value - 2) ** 2)
        assert_allclose(np.mean(values) - 4, -2 * eps + .5 * eps**2, atol=1e-15)


class AnalyzerTests(unittest.TestCase):
    def test_same_corrected_pdf_on_different_fk_grids(self):
        for basis in ('flavor', 'evolution'):
            for enforce in (False, True):
                for hadronic in (False, True):
                    with self.subTest(basis=basis, enforce=enforce, hadronic=hadronic):
                        analyzer = make_analyzer(basis, enforce, hadronic=hadronic)
                        analyzer._evaluate_chi2([0, 1], .15, .1, .2, mode='calibrated')
                        a, b = analyzer.observables
                        assert_allclose(a.last_grid[:, :, -1], b.last_grid[:, :, -1],
                                        rtol=0, atol=1e-14)

    def test_zero_amplitude_has_one_baseline(self):
        for basis in ('flavor', 'evolution'):
            analyzer = make_analyzer(basis)
            baseline = analyzer._evaluate_chi2([], .15, .1, 0, mode='calibrated')
            for coalition in ([0], [1], [0, 1]):
                for sign in (1., -1.):
                    value = analyzer._evaluate_chi2(
                        coalition, .15, .1, 0, mode='calibrated',
                        random_sign_matrix=np.full((1, 2), sign),
                    )
                    assert_allclose(value, baseline, rtol=0, atol=1e-12)

    def test_vector_single_centre_matches_scalar_for_noncontiguous_players(self):
        for basis in ('flavor', 'evolution'):
            for enforce in (False, True):
                scalar = make_analyzer(basis, enforce)
                vector = make_analyzer(basis, enforce, vector=True)
                for coalition in ([1], [0], [0, 1]):
                    for signs in (np.array([[1., -1.]]), np.array([[-1., 1.]])):
                        with self.subTest(basis=basis, enforce=enforce, coalition=coalition):
                            kwargs = dict(mode='calibrated', random_sign_matrix=signs)
                            expected = scalar._evaluate_chi2(coalition, .15, .1, .2, **kwargs)
                            actual = vector._evaluate_chi2(coalition, .15, .1, .2, **kwargs)
                            assert_allclose(actual, expected, rtol=0, atol=1e-12)

    def test_quadratic_shapley_includes_asymmetry_and_cross_terms(self):
        hessian = np.array([[2., .7, -.2], [.7, 3., .4], [-.2, .4, 1.]])
        slope = np.array([1., -2., .5])
        plus, minus = np.array([.1, .2, .3]), np.array([.2, .1, .4])
        coalitions = [c for k in range(4) for c in combinations(range(3), k)]
        cache = {}
        for coalition in coalitions:
            mask = np.zeros(3)
            mask[list(coalition)] = 1
            a, b = mask * plus, -mask * minus
            cache[coalition] = .5 * (slope @ (a + b) + .5 * (a @ hessian @ a + b @ hessian @ b))
        actual = analyzer_module.NNPDFShapleyAnalyzer._compute_shapley_from_cache(cache, coalitions, 3)
        expected = .5 * slope * (plus - minus) + .25 * (
            plus * (hessian @ plus) + minus * (hessian @ minus)
        )
        assert_allclose(actual, expected, atol=1e-15)
        assert_allclose(sum(actual), cache[(0, 1, 2)] - cache[()], atol=1e-15)


if __name__ == '__main__':
    unittest.main()
