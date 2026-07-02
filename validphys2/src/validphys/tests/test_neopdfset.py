"""
Tests for the NeoPDF interpolation backend.
"""

import os

import numpy as np
import pytest

# Skip the whole module if neopdf is not importable
neopdf = pytest.importorskip("neopdf", reason="neopdf not installed")

# In view of making LHAPDF Optional in the future.
try:
    import lhapdf as _lhapdf

    _lhapdf.setVerbosity(0)
    HAS_LHAPDF = True
except ModuleNotFoundError:
    HAS_LHAPDF = False

requires_lhapdf = pytest.mark.skipif(not HAS_LHAPDF, reason="lhapdf not installed")

PDF_NAME = "NNPDF40_nnlo_as_01180"

PIDS = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]

XGRID = np.array([1e-5, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9])
QGRID = np.array([1.7, 5.0, 10.0, 100.0])  # GeV

X_VALUE = 0.1
Q_VALUE = 10.0


@pytest.fixture(scope="module")
def neo_pdfset():
    """LHAPDFSet loaded with the NeoPDF backend."""
    from validphys.lhapdfset import LHAPDFSet

    os.environ["NNPDF_PDF_BACKEND"] = "neopdf"
    pdfset = LHAPDFSet(PDF_NAME, "replicas")
    del os.environ["NNPDF_PDF_BACKEND"]
    return pdfset


@pytest.fixture(scope="module")
def neo_pdfset_t0():
    """LHAPDFSet in t0 mode loaded with the NeoPDF backend."""
    from validphys.lhapdfset import LHAPDFSet

    os.environ["NNPDF_PDF_BACKEND"] = "neopdf"
    pdfset = LHAPDFSet(PDF_NAME, "t0")
    del os.environ["NNPDF_PDF_BACKEND"]
    return pdfset


@pytest.fixture(scope="module")
def lha_pdfset():
    from validphys.lhapdfset import LHAPDFSet

    return LHAPDFSet(PDF_NAME, "replicas")


@pytest.fixture(scope="module")
def lha_pdfset_t0():
    from validphys.lhapdfset import LHAPDFSet

    return LHAPDFSet(PDF_NAME, "t0")


class TestT0Mode:
    """Focused tests on t0 semantics for the NeoPDF backend."""

    def test_t0_members_slice_length(self, neo_pdfset_t0):
        assert len(neo_pdfset_t0.members) == 1

    def test_t0_central_member_identity(self, neo_pdfset_t0):
        assert neo_pdfset_t0.central_member is neo_pdfset_t0.members[0]

    def test_t0_xfxQ_index_zero(self, neo_pdfset_t0):
        val = neo_pdfset_t0.xfxQ(X_VALUE, Q_VALUE, n=0, fl=21)
        assert isinstance(val, float)
        assert np.isfinite(val)

    def test_t0_xfxQ_out_of_range_member_raises(self, neo_pdfset_t0):
        """Requesting member index >= 1 in t0 mode must raise IndexError."""
        with pytest.raises(IndexError):
            neo_pdfset_t0.xfxQ(X_VALUE, Q_VALUE, n=1, fl=21)

    @requires_lhapdf
    def test_t0_xfxQ_agrees_with_lhapdf(self, neo_pdfset_t0, lha_pdfset_t0):
        for fl in [21, 1, 2, -1]:
            if fl not in neo_pdfset_t0.flavors:
                continue
            neo_val = neo_pdfset_t0.xfxQ(X_VALUE, Q_VALUE, n=0, fl=fl)
            lha_val = lha_pdfset_t0.xfxQ(X_VALUE, Q_VALUE, n=0, fl=fl)
            np.testing.assert_equal(neo_val, lha_val)


class TestNeoPDFSetInterface:
    """The following simpy checks that that ``NeoPDFSet`` exposes the same
    interface as ``LHAPDFSet``.
    """

    def test_is_t0_false_for_replicas(self, neo_pdfset):
        assert neo_pdfset.is_t0 is False

    def test_is_t0_true_for_t0(self, neo_pdfset_t0):
        assert neo_pdfset_t0.is_t0 is True

    def test_n_members_positive(self, neo_pdfset):
        assert neo_pdfset.n_members > 0

    def test_n_members_matches_members_length(self, neo_pdfset):
        assert neo_pdfset.n_members == len(neo_pdfset.members)

    def test_t0_has_exactly_one_member(self, neo_pdfset_t0):
        assert neo_pdfset_t0.n_members == 1

    def test_central_member_is_member_zero(self, neo_pdfset):
        assert neo_pdfset.central_member is neo_pdfset.members[0]

    def test_flavors_is_list(self, neo_pdfset):
        assert isinstance(neo_pdfset.flavors, list)

    def test_flavors_contains_gluon(self, neo_pdfset):
        assert 21 in neo_pdfset.flavors

    def test_flavors_contains_quarks(self, neo_pdfset):
        for pid in [-2, -1, 1, 2]:
            assert pid in neo_pdfset.flavors

    def test_flavors_cached(self, neo_pdfset):
        assert neo_pdfset.flavors is neo_pdfset.flavors

    def test_xfxQ_returns_float(self, neo_pdfset):
        val = neo_pdfset.xfxQ(X_VALUE, Q_VALUE, n=0, fl=21)
        assert isinstance(val, float)

    def test_xfxQ_absent_flavour_returns_zero(self, neo_pdfset):
        if 6 not in neo_pdfset.flavors:
            assert neo_pdfset.xfxQ(X_VALUE, Q_VALUE, n=0, fl=6) == 0.0

    def test_grid_values_shape(self, neo_pdfset):
        nx, nq, _ = len(XGRID), len(QGRID), len(PIDS)
        pids = [p for p in PIDS if p in neo_pdfset.flavors]
        result = neo_pdfset.grid_values(np.array(pids), XGRID, QGRID)
        assert result.shape == (neo_pdfset.n_members, len(pids), nx, nq)

    def test_grid_values_dtype_is_float(self, neo_pdfset):
        pids = [p for p in PIDS if p in neo_pdfset.flavors]
        result = neo_pdfset.grid_values(np.array(pids), XGRID, QGRID)
        assert np.issubdtype(result.dtype, np.floating)

    def test_grid_values_finite(self, neo_pdfset):
        pids = [p for p in PIDS if p in neo_pdfset.flavors]
        result = neo_pdfset.grid_values(np.array(pids), XGRID, QGRID)
        assert np.all(np.isfinite(result))

    def test_grid_values_single_x_single_q(self, neo_pdfset):
        pids = [21]
        result = neo_pdfset.grid_values(np.array(pids), np.array([0.1]), np.array([10.0]))
        assert result.shape == (neo_pdfset.n_members, 1, 1, 1)

    def test_t0_grid_values_shape(self, neo_pdfset_t0):
        pids = [p for p in PIDS if p in neo_pdfset_t0.flavors]
        result = neo_pdfset_t0.grid_values(np.array(pids), XGRID, QGRID)
        assert result.shape == (1, len(pids), len(XGRID), len(QGRID))


@requires_lhapdf
class TestNumericalAgreement:
    """
    NeoPDF and LHAPDF must produce identical values (bit-for-bit) for every
    call that both backends support.
    """

    @pytest.mark.parametrize("fl", [21, 1, -1, 2, -2, 3])
    def test_xfxQ_member0(self, neo_pdfset, lha_pdfset, fl):
        if fl not in neo_pdfset.flavors:
            pytest.skip(f"pid {fl} not in set")
        neo_val = neo_pdfset.xfxQ(X_VALUE, Q_VALUE, n=0, fl=fl)
        lha_val = lha_pdfset.xfxQ(X_VALUE, Q_VALUE, n=0, fl=fl)
        np.testing.assert_equal(neo_val, lha_val)

    @pytest.mark.parametrize("fl", [21, 2])
    def test_xfxQ_all_members(self, neo_pdfset, lha_pdfset, fl):
        if fl not in neo_pdfset.flavors:
            pytest.skip(f"pid {fl} not in set")
        for n in range(neo_pdfset.n_members):
            neo_val = neo_pdfset.xfxQ(X_VALUE, Q_VALUE, n=n, fl=fl)
            lha_val = lha_pdfset.xfxQ(X_VALUE, Q_VALUE, n=n, fl=fl)
            np.testing.assert_equal(neo_val, lha_val)

    @pytest.mark.parametrize("x", [1e-5, 1e-3, 0.1, 0.5, 0.9])
    @pytest.mark.parametrize("Q", [1.7, 10.0, 100.0])
    def test_xfxQ_gluon_phase_space(self, neo_pdfset, lha_pdfset, x, Q):
        neo_val = neo_pdfset.xfxQ(x, Q, n=0, fl=21)
        lha_val = lha_pdfset.xfxQ(x, Q, n=0, fl=21)
        np.testing.assert_equal(neo_val, lha_val)

    @pytest.mark.parametrize("x", [1e-5, 0.1, 0.9])
    @pytest.mark.parametrize("Q", [1.7, 100.0])
    @pytest.mark.parametrize("fl", [6, -6])
    def test_xfxQ_absent_flavour_zero(self, neo_pdfset, fl, x, Q):
        """Absent flavours must return 0.0."""
        assert fl not in neo_pdfset.flavors
        assert neo_pdfset.xfxQ(x, Q, n=0, fl=fl) == 0.0

    def test_grid_values_shape_matches_lhapdf(self, neo_pdfset, lha_pdfset):
        pids = np.array([p for p in PIDS if p in neo_pdfset.flavors])
        neo_result = neo_pdfset.grid_values(pids, XGRID, QGRID)
        lha_result = lha_pdfset.grid_values(pids, XGRID, QGRID)
        assert neo_result.shape == lha_result.shape

    @pytest.mark.parametrize("fl", [21, 1, 2])
    def test_grid_values_single_flavour(self, neo_pdfset, lha_pdfset, fl):
        if fl not in neo_pdfset.flavors:
            pytest.skip(f"pid {fl} not in set")
        pids = np.array([fl])
        neo_result = neo_pdfset.grid_values(pids, XGRID, QGRID)
        lha_result = lha_pdfset.grid_values(pids, XGRID, QGRID)
        np.testing.assert_array_equal(neo_result, lha_result)

    def test_grid_values_all_flavours(self, neo_pdfset, lha_pdfset):
        pids = np.array([p for p in PIDS if p in neo_pdfset.flavors])
        neo_result = neo_pdfset.grid_values(pids, XGRID, QGRID)
        lha_result = lha_pdfset.grid_values(pids, XGRID, QGRID)
        np.testing.assert_array_equal(neo_result, lha_result)

    def test_grid_values_member0_equals_xfxQ_scalar(self, neo_pdfset):
        """With nq=1 the grid cell directly maps to the scalar xfxQ value."""
        fl = 21
        x_idx = 2
        x = XGRID[x_idx]
        Q_single = np.array([Q_VALUE])
        grid = neo_pdfset.grid_values(np.array([fl]), XGRID, Q_single)
        scalar = neo_pdfset.xfxQ(x, Q_VALUE, n=0, fl=fl)
        np.testing.assert_equal(grid[0, 0, x_idx, 0], scalar)

    def test_t0_central_value(self, neo_pdfset_t0, lha_pdfset_t0):
        pids = np.array([p for p in PIDS if p in neo_pdfset_t0.flavors])
        neo_result = neo_pdfset_t0.grid_values(pids, XGRID, QGRID)
        lha_result = lha_pdfset_t0.grid_values(pids, XGRID, QGRID)
        np.testing.assert_array_equal(neo_result, lha_result)

    def test_t0_same_as_replica_member0(self, neo_pdfset, neo_pdfset_t0):
        """t0 set must return the same values as member 0 of the replica set."""
        pids = np.array([p for p in PIDS if p in neo_pdfset.flavors])
        full = neo_pdfset.grid_values(pids, XGRID, QGRID)
        t0 = neo_pdfset_t0.grid_values(pids, XGRID, QGRID)
        np.testing.assert_array_equal(full[0:1], t0)

    def test_flavors_match_lhapdf(self, neo_pdfset, lha_pdfset):
        assert sorted(neo_pdfset.flavors) == sorted(lha_pdfset.flavors)

    def test_t0_flavors_match_lhapdf(self, neo_pdfset_t0, lha_pdfset_t0):
        assert sorted(neo_pdfset_t0.flavors) == sorted(lha_pdfset_t0.flavors)

    def test_n_members_matches_lhapdf(self, neo_pdfset, lha_pdfset):
        assert neo_pdfset.n_members == lha_pdfset.n_members


class TestBackendFactory:
    """Verify that ``make_pdf`` in ``lhapdf_compatibility`` selects the correct
    backend based on the ``NNPDF_PDF_BACKEND`` environment variable.
    """

    def test_default_returns_lhapdf_members(self, monkeypatch):
        monkeypatch.delenv("NNPDF_PDF_BACKEND", raising=False)
        from validphys.lhapdf_compatibility import _NeoPDFPDF, make_pdf

        members = make_pdf(PDF_NAME)
        assert not isinstance(members[0], _NeoPDFPDF)

    def test_env_var_neopdf_returns_neopdf_members(self, monkeypatch):
        monkeypatch.setenv("NNPDF_PDF_BACKEND", "neopdf")
        from validphys.lhapdf_compatibility import _NeoPDFPDF, make_pdf

        members = make_pdf(PDF_NAME)
        assert all(isinstance(m, _NeoPDFPDF) for m in members)

    def test_env_var_lhapdf_returns_non_neopdf_members(self, monkeypatch):
        monkeypatch.setenv("NNPDF_PDF_BACKEND", "lhapdf")
        from validphys.lhapdf_compatibility import _NeoPDFPDF, make_pdf

        members = make_pdf(PDF_NAME)
        assert not isinstance(members[0], _NeoPDFPDF)

    def test_invalid_env_var_raises(self, monkeypatch):
        monkeypatch.setenv("NNPDF_PDF_BACKEND", "pdfflow_is_not_valid_here")
        from validphys.lhapdf_compatibility import InvalidPDFBackend, make_pdf

        with pytest.raises(InvalidPDFBackend, match="Unknown backend"):
            make_pdf(PDF_NAME)

    def test_neopdf_single_member_returns_neopdf_member(self, monkeypatch):
        monkeypatch.setenv("NNPDF_PDF_BACKEND", "neopdf")
        from validphys.lhapdf_compatibility import _NeoPDFPDF, make_pdf

        members = make_pdf(PDF_NAME, member=0)
        assert len(members) == 1
        assert isinstance(members[0], _NeoPDFPDF)

    @requires_lhapdf
    def test_both_backends_return_same_member_count(self, monkeypatch):
        from validphys.lhapdf_compatibility import make_pdf

        monkeypatch.setenv("NNPDF_PDF_BACKEND", "neopdf")
        neo_members = make_pdf(PDF_NAME)

        monkeypatch.setenv("NNPDF_PDF_BACKEND", "lhapdf")
        lha_members = make_pdf(PDF_NAME)

        assert len(neo_members) == len(lha_members)
