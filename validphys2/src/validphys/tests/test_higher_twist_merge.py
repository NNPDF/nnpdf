"""Dataset ownership and conditional Gaussian PC posterior transport."""

from types import SimpleNamespace

from matplotlib import pyplot as plt
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import pytest
from scipy.linalg import block_diag

from reportengine.checks import CheckError
from reportengine.resourcebuilder import ResourceError
from validphys.api import API
from validphys.theorycovariance import higher_twist_merge as merge
from validphys.theorycovariance import higher_twist_posterior as pc


def _fit(name, dmean, dcov, regression, residual, jmean, factorised=False):
    """Generate a joint Gaussian independently via J = muJ + A(D-muD) + eps."""
    parameters = {"f2p": {"nodes": [0, 0.5, 1], "yshift": [1.0, 1.0, 0.0]}}
    for key in pc._DIJET_FACTORS if factorised else ("Hj",):
        parameters[key] = {"nodes": [0, 1], "yshift": [1.0, 1.0]}
    basis = pc.pc_parameter_basis(parameters)
    nj = len(jmean)
    transform = np.eye(3 + nj)
    transform[3:, :2] = regression
    covariance = transform @ block_diag(dcov, np.ones((1, 1)), residual) @ transform.T
    posterior = pc.PCPosterior(
        pd.Series(np.r_[dmean, 0.0, jmean], index=basis.columns),
        pd.DataFrame(covariance, index=basis.columns, columns=basis.columns),
        pd.DataFrame(np.eye(len(basis.columns)), index=basis.columns, columns=basis.columns),
    )
    # Shared DIS, a private jet sector, and an unaffected dataset. Nonconsecutive IDs.
    index = pd.MultiIndex.from_tuples(
        [
            ("DIS NC", "DIS", 2),
            ("DIS NC", "DIS", 7),
            ("JETS", name, 5),
            ("JETS", name, 9),
            ("DY", "unchanged", 3),
        ],
        names=["group", "dataset", "id"],
    )
    beta = np.zeros((5, len(basis.columns)))
    beta[:2, :2] = [[1, 0.2], [0.3, 1.1]]
    beta[2:4, 3:] = np.arange(1, 2 * nj + 1).reshape(2, nj) / 5
    # Node response used only for a simple physical mean in the merge tests.
    node_response = np.zeros((5, len(basis)))
    node_response[:2, :2] = beta[:2, :2]
    node_response[2:4, 3:] = 0.5
    return merge.PCMergeFit(
        name,
        parameters,
        basis,
        posterior,
        pd.DataFrame(beta, index=index, columns=basis.columns),
        node_response,
        None,
    )


@pytest.fixture
def inputs(monkeypatch):
    items = [
        _fit(
            "reference",
            [1, -1],
            [[2, 0.4], [0.4, 1]],
            [[0.2, -0.3], [0.4, 0.6]],
            [[1, 0.2], [0.2, 2]],
            [2, 3],
        ),
        _fit(
            "sector2",
            [-1, 2],
            [[1, -0.1], [-0.1, 3]],
            [[-0.5, 0.1], [0.2, 0.7]],
            [[2, -0.3], [-0.3, 1]],
            [-2, 1],
        ),
        _fit(
            "sector3",
            [0.5, 0.2],
            [[3, 0.2], [0.2, 2]],
            [[0.3, 0.2], [-0.1, 0.4]],
            [[1, 0.1], [0.1, 1]],
            [1, -3],
        ),
    ]

    def mean_shift(datasets, pdf, parameters, posterior, index):
        full_index = next(item.shifts.index for item in items if item.datasets is datasets)
        return pd.Series(datasets @ posterior.mean, index=full_index, name="mean_shift").loc[index]

    monkeypatch.setattr(pc, "_prediction_mean", mean_shift)
    return items


@pytest.mark.parametrize("factor", [0, 1, 2, 3.5])
def test_independent_merge_preserves_owned_blocks(inputs, factor):
    result = API.pc_merged_posterior(pc_merge_inputs=inputs)
    indices = [inputs[0].shifts.index, inputs[1].shifts.index[2:4], inputs[2].shifts.index[2:4]]
    expected_cov, expected_mean = [], []
    for item, index in zip(inputs, indices):
        expected_cov.append(pc.pc_prediction_covmat(item.posterior, item.shifts).loc[index, index])
        expected_mean.append(
            pd.Series(
                item.datasets @ pc.pc_posterior(item.posterior, item.basis).mean,
                index=item.shifts.index,
            ).loc[index]
        )
    mean = np.concatenate(expected_mean)
    covariance = block_diag(*expected_cov)
    assert_allclose(result.mean, mean)
    assert_allclose(result.covariance, covariance)
    output = API.pc_merged_covmat(pc_merged_posterior=result, pc_covariance_factor=factor)
    assert_allclose(output, factor * covariance + np.outer(mean, mean))
    assert_allclose(output.iloc[:5, 5:], np.outer(mean[:5], mean[5:]))
    assert list(result.owners) == ["reference"] * 5 + ["sector2"] * 2 + ["sector3"] * 2
    assert result.mean.index.droplevel(0).is_unique
    assert np.linalg.eigvalsh(output).min() > -1e-10


def test_conditional_merge_matches_explicit_joint_gaussian(inputs):
    result = API.pc_merged_posterior(pc_merge_inputs=inputs, pc_merge_method="conditional")
    dcov = inputs[0].posterior.covariance.to_numpy()[:2, :2]
    dmean = inputs[0].posterior.mean.to_numpy()[:2]
    regressions, residuals, means, responses = [], [], [], []
    for n, item in enumerate(inputs):
        covariance = item.posterior.covariance.to_numpy()
        regression = covariance[3:, :2] @ np.linalg.inv(covariance[:2, :2])
        regressions.append(regression)
        residuals.append(covariance[3:, 3:] - regression @ covariance[:2, 3:])
        means.append(
            item.posterior.mean.to_numpy()[3:]
            + regression @ (dmean - item.posterior.mean.to_numpy()[:2])
        )
        responses.append(
            item.shifts.to_numpy()[:, [0, 1, 3, 4]] if n == 0 else item.shifts.to_numpy()[2:4, 3:]
        )
    transform = np.eye(8)
    transform[2:, :2] = np.vstack(regressions)
    joint = transform @ block_diag(dcov, *residuals) @ transform.T
    response = block_diag(*responses)
    assert_allclose(result.covariance, response @ joint @ response.T, atol=1e-13)
    # Reference block is unchanged, including its DIS--jet correlations.
    original = pc.pc_prediction_covmat(inputs[0].posterior, inputs[0].shifts)
    assert_allclose(result.covariance.iloc[:5, :5], original, atol=1e-13)
    # Means are physical node-response projections of the transported means.
    node_responses = [inputs[0].datasets[:, [0, 1, 3, 4]]] + [
        i.datasets[2:4, 3:] for i in inputs[1:]
    ]
    assert_allclose(result.mean, block_diag(*node_responses) @ np.concatenate([dmean, *means]))
    assert np.linalg.eigvalsh(result.covariance).min() > -1e-10
    assert np.any(result.covariance.iloc[:2, 5:].to_numpy() != 0)


def test_factorised_transport_keeps_full_combination_covariance():
    item = _fit(
        "dijets",
        [1, 2],
        [[2, 0.1], [0.1, 1]],
        [[0.1, 0.2], [0.2, -0.3], [0.3, 0.4], [-0.4, 0.2]],
        np.eye(4),
        [1, 2, -1, 0],
        factorised=True,
    )
    reference = _fit(
        "reference", [-1, 1], [[1, 0.2], [0.2, 3]], np.zeros((2, 2)), np.eye(2), [0, 0]
    )
    result, residual, gain, dd = merge._conditional_transport(item, reference)
    assert np.linalg.matrix_rank(item.basis) < item.basis.shape[1]
    # Transporting combinations and then projecting gives the direct node conditional.
    nodes = pc.pc_posterior(item.posterior, item.basis)
    cov = nodes.covariance.to_numpy()
    regression = cov[:, :2] @ np.linalg.inv(cov[:2, :2])
    expected = cov - regression @ cov[:2, :] + regression @ dd @ regression.T
    assert_allclose(pc.pc_posterior(result, item.basis).covariance, expected, atol=1e-13)
    assert_allclose(
        item.basis.to_numpy() @ result.mean,
        nodes.mean
        + regression
        @ (reference.posterior.mean.iloc[:2].to_numpy() - item.posterior.mean.iloc[:2].to_numpy()),
        atol=1e-13,
    )
    same, _, _, _ = merge._conditional_transport(item, item)
    assert_allclose(
        pc.pc_prediction_covmat(same, item.shifts),
        pc.pc_prediction_covmat(item.posterior, item.shifts),
        atol=1e-13,
    )
    assert np.linalg.eigvalsh(residual).min() > -1e-12


def test_conditional_equals_independent_when_dis_cross_covariances_vanish(inputs):
    items = []
    for item in inputs:
        covariance = item.posterior.covariance.copy()
        covariance.iloc[:3, 3:] = 0
        covariance.iloc[3:, :3] = 0
        items.append(item._replace(posterior=item.posterior._replace(covariance=covariance)))
    independent = merge.pc_merged_posterior(items)
    conditional = merge.pc_merged_posterior(items, "conditional")
    assert_allclose(conditional.mean, independent.mean)
    assert_allclose(conditional.covariance, independent.covariance, atol=1e-13)


def test_ownership_reference_selection_and_subset(inputs):
    result = merge.pc_merged_posterior(
        inputs, pc_merge_reference="sector2", pc_merge_datasets=["DIS", "sector3"]
    )
    assert list(result.owners) == ["sector2", "sector2", "sector3", "sector3"]
    assert list(result.mean.index.get_level_values("id")) == [2, 7, 5, 9]
    result = merge.pc_merged_posterior(
        inputs, pc_dataset_owners={"DIS": "sector3"}, pc_merge_datasets=["DIS"]
    )
    assert set(result.owners) == {"sector3"}
    for kwargs, match in (
        ({"pc_dataset_owners": {"DIS": "missing"}}, "not available"),
        ({"pc_merge_reference": "missing"}, "reference fit"),
        ({"pc_merge_datasets": ["unknown"]}, "unavailable"),
        ({"pc_merge_datasets": []}, "No owned"),
        ({"pc_merge_method": "conditional", "pc_dataset_owners": {"DIS": "sector2"}}, "all DIS"),
    ):
        with pytest.raises(CheckError, match=match):
            merge.pc_merged_posterior(inputs, **kwargs)


def test_dis_alignment_and_fixed_nodes(inputs):
    reference, item = inputs[:2]
    reordered = item._replace(basis=item.basis.iloc[::-1])
    expected = merge._conditional_transport(item, reference)[0]
    actual = merge._conditional_transport(reordered, reference)[0]
    assert_allclose(actual.mean, expected.mean)
    assert_allclose(actual.covariance, expected.covariance)
    for par, match in (
        ({"nodes": [0, 0.4, 1], "yshift": [1, 1, 0]}, "positions differ"),
        ({"nodes": [0, 0.5, 1], "yshift": [1, 1, 1]}, "Fixed DIS nodes differ"),
    ):
        changed = item._replace(parameters={**item.parameters, "f2p": par})
        with pytest.raises(CheckError, match=match):
            merge._conditional_transport(changed, reference)
    comparison = merge.pc_merge_dis_comparison(inputs)
    assert_allclose(comparison.loc["reference", "difference_over_combined_std"], 0)
    assert_allclose(comparison.loc[("sector2", "f2p", "f2p(2)"), "std"], 0)
    matrices = merge.pc_merge_dis_covmats(inputs)
    assert_allclose(
        matrices.loc["reference"],
        reference.posterior.covariance.iloc[:3, :3] * np.outer([1, 1, 0], [1, 1, 0]),
    )


@pytest.mark.parametrize("factor", [-1, np.inf, np.nan, "two", True])
def test_invalid_factor_is_rejected_before_execution(factor):
    with pytest.raises(ResourceError, match="non-negative"):
        API.pc_merged_covmat(pc_merged_posterior=None, pc_covariance_factor=factor)


def test_invalid_merge_method_is_rejected():
    with pytest.raises(ResourceError, match="pc_merge_method"):
        API.pc_merged_posterior(pc_merge_inputs=[], pc_merge_method="other")


def test_nested_fit_resolution_preserves_environment():
    # In particular, a CLI --no-net loader must not become a default network loader.
    assert API.pc_merge_api().loadedenv is API.loadedenv


@pytest.mark.parametrize("method", ["independent", "conditional"])
def test_merged_heatmaps_show_dataset_labels_in_export_order(inputs, method):
    posterior = merge.pc_merged_posterior(inputs, pc_merge_method=method)
    matrix = merge.pc_merged_covmat(posterior, pc_covariance_factor=3)
    active = np.diag(matrix) > 0
    expected = matrix.loc[active, active]
    covfig = API.plot_pc_merged_covmat(pc_merged_covmat=matrix)
    corrfig = API.plot_pc_merged_corrmat(pc_merged_covmat=matrix)
    try:
        assert_allclose(covfig.axes[0].images[0].get_array(), expected)
        assert_allclose(
            corrfig.axes[0].images[0].get_array(),
            expected / np.sqrt(np.outer(np.diag(expected), np.diag(expected))),
        )
        assert [label.get_text() for label in covfig.axes[0].get_xticklabels()] == [
            "DIS",
            "reference",
            "sector2",
            "sector3",
        ]
    finally:
        plt.close(covfig)
        plt.close(corrfig)
    # A valid zero error model must still render, rather than fail on empty axes/log bounds.
    zero = matrix * 0
    fig = merge.plot_pc_merged_covmat(zero)
    try:
        fig.canvas.draw()
        assert_allclose(fig.axes[0].images[0].get_array(), zero)
    finally:
        plt.close(fig)


def test_prediction_mean_product_and_ratio_conventions(monkeypatch):
    parameters = {key: {"nodes": [0, 1], "yshift": [1, 1]} for key in pc._DIJET_FACTORS}
    basis = pc.pc_parameter_basis(parameters)
    posterior = pc.PCPosterior(
        pd.Series([1, 2, 3, 4], index=basis.index),
        pd.DataFrame(np.eye(4) + 0.2, index=basis.index, columns=basis.index),
        basis @ basis.T,
    )
    index = pd.MultiIndex.from_tuples(
        [("DIJET", "jets", 3), ("OTHER", "zero", 0)], names=["group", "dataset", "id"]
    )

    def response(ds, pdf, pars, combinations=pc.construct_pars_combs):
        return {
            c["label"]: np.array(
                [c["comb"][pc._DIJET_FACTORS[0]][0] * c["comb"][pc._DIJET_FACTORS[1]][1]]
            )
            for c in combinations(pars, pc._DIJET_FACTORS)
        }

    monkeypatch.setattr(pc, "compute_deltas_pc", response)
    result = pc.pc_prediction_mean(
        [SimpleNamespace(name="jets")], None, parameters, posterior, index
    )
    assert_allclose(result, [1 * 4 + 0.2, 0])
    # For a ratio, the existing convention is evaluation at means, not B @ mean.
    monkeypatch.setattr(
        pc,
        "compute_deltas_pc",
        lambda ds, pdf, pars, comb: {
            "all": np.array([(1 + pars["f2d"]["yshift"][0]) / (1 + pars["f2p"]["yshift"][0]) - 1])
        },
    )
    parameters = {key: {"nodes": [0, 1], "yshift": [1, 1]} for key in ("f2p", "f2d")}
    basis = pc.pc_parameter_basis(parameters)
    posterior = posterior._replace(mean=pd.Series([1, 2, 3, 4], index=basis.index))
    result = pc.pc_prediction_mean(
        [SimpleNamespace(name="jets")], None, parameters, posterior, index
    )
    assert_allclose(result, [1, 0])
