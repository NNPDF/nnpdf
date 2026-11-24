"""Numerical conditioning, data alignment and FK-level PC export."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import pytest

from reportengine.checks import CheckError
from validphys.api import API
from validphys.fkparser import parse_cfactor
from validphys.theorycovariance import higher_twist_posterior as pc


def test_condition_matches_parameter_space_update():
    """Independent precision-matrix calculation, with an exactly fixed node."""
    rng = np.random.default_rng(51)
    response = rng.normal(size=(7, 3))
    basis = np.diag([1.2, 0.8, 0.0])
    beta = response @ basis
    noise = rng.normal(size=(7, 7))
    experimental = noise @ noise.T + np.eye(7)
    predictions = rng.normal(size=(7, 19))
    data = rng.normal(size=7)
    mean, cov, prior = pc._condition(beta, basis, experimental + beta @ beta.T, predictions, data)
    response = response[:, :2]
    inv_exp = np.linalg.inv(experimental)
    conditional = np.linalg.inv(np.diag([1 / 1.2**2, 1 / 0.8**2]) + response.T @ inv_exp @ response)
    gain = conditional @ response.T @ inv_exp
    assert_allclose(mean[:2], gain @ (data - predictions.mean(axis=1)))
    expected = conditional + gain @ np.cov(predictions, bias=True) @ gain.T
    assert_allclose(cov[:2, :2], expected)
    assert_allclose(cov[2], 0, atol=1e-15)
    assert mean[2] == 0
    assert_allclose(prior, basis @ basis.T)


def test_factorised_parameter_basis():
    parameters = {
        "H2j_ystar": {"nodes": [0, 1], "yshift": [2.0, 3.0]},
        "H2j_yb": {"nodes": [0, 1], "yshift": [4.0, 5.0]},
    }
    basis = pc.pc_parameter_basis(parameters)
    assert_allclose(basis, [[2, 2, 0, 0], [0, 0, 3, 3], [4, 0, 4, 0], [0, 5, 0, 5]])
    # The factorised convention has a singular prior, which conditioning must allow.
    assert np.linalg.matrix_rank(basis) == 3
    beta = np.array([[1, 2, 3, 4], [-1, 0, 1, 2]])
    _, covariance, _ = pc._condition(
        beta, basis.to_numpy(), np.eye(2) + beta @ beta.T, np.zeros((2, 3)), np.ones(2)
    )
    assert np.linalg.eigvalsh(covariance).min() > -1e-12


def test_condition_is_invariant_under_data_units():
    """DIS and hadronic cross sections span many orders of magnitude."""
    rng = np.random.default_rng(71)
    beta = rng.normal(size=(50, 3))
    noise = rng.normal(size=(50, 10))
    covariance = np.eye(50) + noise @ noise.T + beta @ beta.T
    predictions = rng.normal(size=(50, 12))
    data = rng.normal(size=50)
    reference = pc._condition(beta, np.eye(3), covariance, predictions, data)
    units = np.geomspace(1e-9, 1e9, 50)
    actual = pc._condition(
        units[:, None] * beta,
        np.eye(3),
        units[:, None] * covariance * units,
        units[:, None] * predictions,
        units * data,
    )
    for result, expected in zip(actual, reference):
        assert_allclose(result, expected, rtol=1e-11, atol=1e-12)


@pytest.mark.parametrize("diagonal", [False, True])
def test_pseudodata_alignment(diagonal):
    """Both group order and training/validation ordering differ from predictions."""
    rng = np.random.default_rng(10)
    n = 12  # exercise numeric rather than lexicographic eigenmode ordering
    data_index = pd.MultiIndex.from_tuples(
        [("experiment", "A" if i < 6 else "B", i) for i in range(n)],
        names=["experiment", "dataset", "id"],
    )
    groups_index = pd.MultiIndex.from_tuples(
        [("process", ds, i) for _, ds, i in reversed(data_index)], names=["group", "dataset", "id"]
    )
    samples = rng.normal(size=(3, n))
    rotation = rng.normal(size=(n, n)) + 10 * np.eye(n) if diagonal else None
    index = [f"eigenmode {i}" for i in range(n)] if diagonal else data_index
    frames = []
    for sample in samples:
        values = rotation @ sample if diagonal else sample
        frame = pd.DataFrame(values, index=index).iloc[rng.permutation(n)]
        frames.append(SimpleNamespace(pseudodata=frame))
    result = pc.pc_central_pseudodata(frames, rotation, data_index, groups_index)
    assert_allclose(result, samples.mean(axis=0)[::-1])
    assert result.index.equals(groups_index)


def test_covmat_alignment(tmp_path):
    index = pd.MultiIndex.from_tuples(
        [("old", "A", 2), ("old", "B", 1)], names=["group", "dataset", "id"]
    )
    matrix = pd.DataFrame([[2.0, 0.4], [0.4, 3.0]], index=index, columns=index)
    path = tmp_path / "covmat.csv"
    matrix.iloc[::-1, :].to_csv(path, sep="\t")
    target = pd.MultiIndex.from_tuples([("new", "B", 1), ("new", "A", 2)], names=index.names)
    assert_allclose(pc._read_covmat(path, target), [[3, 0.4], [0.4, 2]])
    with pytest.raises(KeyError):
        pc._read_covmat(path, pd.MultiIndex.from_tuples([("new", "missing", 0)]))


def test_posterior_uses_other_theory_uncertainties(tmp_path):
    index = pd.MultiIndex.from_tuples([("DIS", "test", 0)], names=["group", "dataset", "id"])
    parameters = {"f2p": {"nodes": [0, 1], "yshift": [2.0, 0.0]}}
    basis = pc.pc_parameter_basis(parameters)
    beta = pd.DataFrame([[2.0, 0.0]], index=index, columns=basis.columns)
    tables = tmp_path / "tables"
    tables.mkdir()
    prefix = "datacuts_theory_theorycovmatconfig_"
    pc_path = tables / f"{prefix}point_prescriptions1_theory_covmat_custom_per_prescription.csv"
    for path, value in (
        (tables / f"{prefix}point_prescriptions0_theory_covmat_custom_per_prescription.csv", 3),
        (pc_path, 4),
        (tables / f"{prefix}user_covmat.csv", 2),
    ):
        pd.DataFrame([[value]], index=index, columns=index).to_csv(path, sep="\t")
    fit = SimpleNamespace(
        path=tmp_path,
        as_input=lambda: {
            "theorycovmatconfig": {
                "point_prescriptions": ["7 point", "power corrections"],
                "user_covmat_path": "user.csv",
            }
        },
    )
    prediction_table = pd.DataFrame(
        [[99, 99, 1, 3]],
        index=index,
        columns=["data_central", "theory_central", "rep_00001", "rep_00002"],
    )
    experimental = pd.DataFrame([[1.0]], index=index, columns=index)
    central_data = pd.Series([0.0], index=index)
    combinations = pc.pc_combination_posterior(
        fit, beta, experimental, prediction_table, central_data, [7, 2]
    )
    posterior = pc.pc_posterior(combinations, basis)
    # Total variance = 1 experimental + 3 MHOU + 4 PC + 2 user; PDF variance = 1.
    assert_allclose(posterior.mean, [-0.8, 0])
    assert_allclose(posterior.covariance, [[4 - 16 / 10 + 16 / 100, 0], [0, 0]])
    # Exercise the registered table provider through reportengine too.
    result = API.pc_posterior_table(pc_posterior=posterior, pc_parameters=parameters)
    assert_allclose(result["mean"], posterior.mean)
    pd.DataFrame([[8]], index=index, columns=index).to_csv(pc_path, sep="\t")
    with pytest.raises(CheckError, match="disagree"):
        pc.pc_combination_posterior(fit, beta, experimental, prediction_table, central_data, [7, 2])


def test_posterior_preserves_cross_process_experimental_correlations(tmp_path):
    index = pd.MultiIndex.from_tuples(
        [("DIS NC", "A", 0), ("DIS CC", "B", 0)], names=["group", "dataset", "id"]
    )
    parameters = {"f2p": {"nodes": [0, 1], "yshift": [1.0, 0.0]}}
    basis = pc.pc_parameter_basis(parameters)
    beta = pd.DataFrame([[1.0, 0.0], [0.0, 0.0]], index=index, columns=basis.columns)
    tables = tmp_path / "tables"
    tables.mkdir()
    path = tables / (
        "datacuts_theory_theorycovmatconfig_point_prescriptions0"
        "_theory_covmat_custom_per_prescription.csv"
    )
    pd.DataFrame([[1.0, 0.0], [0.0, 0.0]], index=index, columns=index).to_csv(path, sep="\t")
    fit = SimpleNamespace(
        path=tmp_path,
        as_input=lambda: {"theorycovmatconfig": {"point_prescriptions": ["power corrections"]}},
    )
    # Both process groups belong to the same experiment; the covariance arrives
    # in reverse dataset order, with a different top-level label.
    exp_index = pd.MultiIndex.from_tuples([("experiment", "B", 0), ("experiment", "A", 0)])
    experimental = pd.DataFrame([[2.0, 0.5], [0.5, 1.0]], index=exp_index, columns=exp_index)
    predictions = pd.DataFrame(0.0, index=index, columns=["rep_00001", "rep_00002"])
    data = pd.Series([0.0, 1.0], index=index)
    combinations = pc.pc_combination_posterior(fit, beta, experimental, predictions, data, [1, 2])
    result = pc.pc_posterior(combinations, basis)
    # Only B has a residual. Its correlation with A must still shift A's PC.
    assert_allclose(result.mean, [-0.5 / 3.75, 0.0])
    assert_allclose(result.covariance, [[1 - 2 / 3.75, 0.0], [0.0, 0.0]])


@pytest.mark.parametrize("factorised", [False, True])
def test_prediction_space_projection_matches_data_space_update(factorised):
    """Compare B P_combination B.T to the independent full data-space formula."""
    rng = np.random.default_rng(18)
    keys = pc._DIJET_FACTORS if factorised else ("Hj",)
    parameters = {key: {"nodes": [0, 1], "yshift": [2.0, 3.0]} for key in keys}
    basis = pc.pc_parameter_basis(parameters)
    index = pd.MultiIndex.from_tuples(
        [("JETS", "A", i) for i in range(3)] + [("JETS", "B", i) for i in range(4)],
        names=["group", "dataset", "id"],
    )
    beta = rng.normal(size=(7, basis.shape[1]))
    beta[-1] = 0  # an unaffected point retains a zero covariance row
    noise = rng.normal(size=(7, 7))
    s = beta @ beta.T
    total = noise @ noise.T + np.eye(7) + s
    predictions = rng.normal(size=(7, 12))
    predictions[0] = 0  # relative uncertainty is undefined here
    data = rng.normal(size=7)
    mean, cov, prior = pc._condition(beta, np.eye(beta.shape[1]), total, predictions, data)
    combinations = pc.PCPosterior(
        pd.Series(mean, index=basis.columns),
        pd.DataFrame(cov, index=basis.columns, columns=basis.columns),
        pd.DataFrame(prior, index=basis.columns, columns=basis.columns),
    )
    shifts = pd.DataFrame(beta, index=index, columns=basis.columns)
    result = pc.pc_prediction_covmat(combinations, shifts)
    inv = np.linalg.inv(total)
    expected = s - s @ inv @ s + s @ inv @ np.cov(predictions, bias=True) @ inv @ s
    assert_allclose(result, expected, rtol=1e-12, atol=1e-12)
    assert result.index.equals(index) and result.columns.equals(index)
    assert_allclose(result.iloc[-1], 0)
    correlations = pc.pc_prediction_corrmat(result)
    assert_allclose(np.diag(correlations), [1, 1, 1, 1, 1, 1, 0])
    assert_allclose(correlations.iloc[-1], 0)
    # Sharing the combination posterior must preserve the existing node update.
    nodes = pc.pc_posterior(combinations, basis)
    direct = pc._condition(beta, basis.to_numpy(), total, predictions, data)
    assert_allclose(nodes.mean, direct[0], atol=1e-12)
    assert_allclose(nodes.covariance, direct[1], atol=1e-12)
    table = pd.DataFrame(predictions, index=index, columns=[f"rep_{i:05}" for i in range(12)])
    uncertainties = pc.pc_prediction_uncertainties(combinations, shifts, table.iloc[::-1])
    assert_allclose(uncertainties["posterior_std"] ** 2, np.diag(expected), atol=1e-12)
    assert_allclose(uncertainties["prior_std"] ** 2, np.diag(s))
    assert np.isnan(uncertainties["posterior_relative_std"].iloc[0])
    assert_allclose(
        uncertainties["posterior_relative_std"].iloc[1:],
        np.sqrt(np.maximum(0, np.diag(expected)[1:])) / abs(predictions.mean(axis=1)[1:]),
    )


def test_correlated_interpolation_and_product():
    weights = np.array([[0.5, 0.5, 0, 0]])
    average = np.array([0, 0, 0.5, 0.5])
    mean = np.array([1, 3, 2, 4])
    covariance = np.eye(4) + np.ones((4, 4)) * 0.2
    central, std = pc._curve_moments(weights, mean, covariance)
    assert_allclose(central, [2])
    assert_allclose(std**2, [0.7])
    central, std = pc._curve_moments(weights, mean, covariance, average)
    assert_allclose(central, [6.2])
    assert_allclose(std**2, [4 * 0.7 + 9 * 0.7 + 2 * 2 * 3 * 0.2 + 0.7**2 + 0.2**2])


def test_experimental_comparison_aligns_points_and_uses_data_normalisation():
    """Experiment order differs; data may be negative/zero and errors may be zero."""
    index = pd.MultiIndex.from_tuples(
        [("JETS", "A", 7), ("DIS NC", "B", 2), ("JETS", "A", 9), ("DY", "C", 0)],
        names=["group", "dataset", "id"],
    )
    points = pd.DataFrame({"prior_std": [4, 5, 6, 0], "posterior_std": [2, 3, 4, 0]}, index=index)
    exp_index = pd.MultiIndex.from_tuples(
        [("experiment", ds, i) for _, ds, i in index], names=index.names
    )
    covariance = pd.DataFrame(np.diag([16.0, 36, 0, 1]), index=exp_index, columns=exp_index)
    central = pd.DataFrame(
        {"data_central": [-10.0, 0, 20, 1], "theory_central": [99.0] * 4}, index=index
    )
    result = API.pc_prediction_experimental_uncertainties(
        pc_prediction_uncertainties=points,
        experiments_covmat_no_table=covariance.iloc[::-1, ::-1],
        group_result_table_no_table=central.iloc[::-1],
    )
    assert result.index.equals(index)
    assert_allclose(result["experimental_std"], [4, 6, 0, 1])
    assert_allclose(result["posterior_relative_std"], [0.2, np.nan, 0.2, 0])
    assert_allclose(result["experimental_relative_std"], [0.4, np.nan, 0, 1])
    assert_allclose(result["posterior_over_experimental"], [0.5, 0.5, np.nan, 0])


def _dataset(name, process, kinematics, specs=()):
    return SimpleNamespace(
        name=name,
        fkspecs=specs,
        cuts=SimpleNamespace(load=lambda: np.array([0])),
        commondata=SimpleNamespace(
            metadata=SimpleNamespace(
                process_type=SimpleNamespace(name=process),
                experiment="CMS",
                load_kinematics=lambda: pd.DataFrame(kinematics),
            )
        ),
    )


def _posterior(parameters, mean, covariance=None):
    basis = pc.pc_parameter_basis(parameters)
    index = basis.index
    if covariance is None:
        covariance = np.eye(len(index))
    return pc.PCPosterior(
        pd.Series(mean, index=index),
        pd.DataFrame(covariance, index=index, columns=index),
        basis @ basis.T,
    )


def test_ratio_cfactor_operands_without_pdf():
    parameters = {key: {"nodes": [0, 1], "yshift": [1, 1]} for key in ("f2p", "f2d")}
    posterior = _posterior(parameters, [2, 4, 6, 10])
    ds = _dataset("NMC_NC_NOTFIXED_EM-F2", "DIS_NC", {"x": [0.25, 0.5], "Q2": [5, 10]})
    numerator, denominator = pc._posterior_cfactor(ds, parameters, posterior)
    assert_allclose(numerator, 1 + np.array([7, 8]) / [5, 10])
    assert_allclose(denominator, 1 + np.array([2.5, 3]) / [5, 10])
    assert ds.cuts is not None


@pytest.mark.parametrize(
    "process,kinematics",
    [
        ("DIJET_3D", {"ystar": [0.25], "yb": [0.6], "m_jj": [100]}),
        ("DIJET", {"ymax": [0.8], "m_jj": [100]}),
    ],
)
def test_factorised_cfactor_includes_covariance(process, kinematics):
    parameters = {key: {"nodes": [0, 1], "yshift": [1, 1]} for key in pc._DIJET_FACTORS}
    covariance = np.eye(4) + np.ones((4, 4)) * 0.2
    ds = _dataset("test", process, kinematics)
    # With zero means the entire correction comes from the cross-covariance.
    posterior = _posterior(parameters, [0, 0, 0, 0], covariance)
    (factor,) = pc._posterior_cfactor(ds, parameters, posterior)
    integration_length = 1 if process == "DIJET_3D" else 0.8
    assert_allclose(factor, [1 + 0.2 * integration_length / 100])


def test_cfactor_export_split_tables_and_missing_bins(tmp_path, monkeypatch):
    import pineappl.fk_table

    paths = [Path("a.pineappl.lz4"), Path("b.pineappl.lz4")]
    sizes = dict(zip(paths, [2, 3]))
    monkeypatch.setattr(
        pineappl.fk_table,
        "FkTable",
        SimpleNamespace(
            read=lambda path: SimpleNamespace(bin_normalizations=lambda: np.ones(sizes[path]))
        ),
    )
    spec = SimpleNamespace(
        fkpath=paths, legacy=False, metadata=SimpleNamespace(shifts={"a": 1, "b": 1})
    )
    ds = _dataset("test", "JET", {"y": np.linspace(0, 1, 7), "pT": np.arange(7) + 10}, [spec])
    parameters = {"Hj": {"nodes": [0, 1], "yshift": [1, 1]}}
    posterior = _posterior(parameters, [2, 2])
    manifest = pc.pc_posterior_cfactors(
        SimpleNamespace(name="testfit"), [ds], parameters, posterior, tmp_path
    )
    for row, indices in zip(manifest.itertuples(), [[1, 2], [4, 5, 6]]):
        with (tmp_path / row.path).open("rb") as stream:
            result = parse_cfactor(stream)
        assert_allclose(result.central_value, 1 + 2 / (np.array(indices) + 10))
        assert_allclose(result.uncertainty, 0)
    assert ds.cuts is not None


def test_fixed_nodes_have_zero_correlations():
    parameters = {"Hj": {"nodes": [0, 1], "yshift": [1, 0]}}
    posterior = _posterior(parameters, [2, 0], np.diag([2.0, 0]))
    assert_allclose(pc.pc_posterior_corrmat(posterior), [[1, 0], [0, 0]])
