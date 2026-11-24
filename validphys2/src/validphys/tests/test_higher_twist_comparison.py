"""Central jet-correction criterion, fitted point alignment and combination-space bands."""

from types import SimpleNamespace

from matplotlib import pyplot as plt
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import pytest

from reportengine.checks import CheckError
from validphys.api import API
from validphys.theorycovariance import higher_twist_comparison as comparison
from validphys.theorycovariance import higher_twist_posterior as pc


def _inputs(process, kinematics, parameters, mean, covariance, central, ids):
    basis = pc.pc_parameter_basis(parameters)
    combinations = pc.PCPosterior(
        pd.Series(mean, index=basis.columns),
        pd.DataFrame(covariance, index=basis.columns, columns=basis.columns),
        pd.DataFrame(np.eye(len(mean)), index=basis.columns, columns=basis.columns),
    )
    commondata = SimpleNamespace(central_values=pd.Series(central))
    dataset = SimpleNamespace(
        name="test",
        cuts=SimpleNamespace(load=lambda: np.sort(ids)),
        commondata=SimpleNamespace(
            load=lambda: commondata,
            metadata=SimpleNamespace(
                process_type=SimpleNamespace(name=process),
                experiment="CMS",
                load_kinematics=lambda: pd.DataFrame(kinematics),
            ),
        ),
    )
    index = pd.MultiIndex.from_tuples(
        [("JETS", "test", i) for i in ids], names=["group", "dataset", "id"]
    )
    exp_index = pd.MultiIndex.from_tuples([("CMS", "test", i) for i in ids], names=index.names)
    experimental = pd.DataFrame(np.eye(len(ids)), index=exp_index, columns=exp_index)
    return dict(
        pc_datasets=[dataset],
        pc_parameters=parameters,
        pc_posterior=pc.pc_posterior(combinations, basis),
        pc_combination_posterior=combinations,
        experiments_covmat_no_table=experimental,
        groups_index=index,
    )


@pytest.fixture
def jets():
    args = _inputs(
        "JET",
        {"y": [0, 0.5, 1, 0, 0.5, 1], "pT": [10, 20, 40, 80, 160, 320]},
        {"Hj": {"nodes": [0, 1], "yshift": [2, 4]}},
        [-2, -1],
        [[0.25, 0.1], [0.1, 0.5]],
        [-10, -20, 99, 0, 20, 10],
        [5, 1, 4, 3],
    )
    exp = args["experiments_covmat_no_table"]
    exp.iloc[:, :] = np.diag(np.array([0.05, 4, 0.8, 2]) ** 2)
    args["experiments_covmat_no_table"] = exp.iloc[::-1, ::-1]
    combinations = args["pc_combination_posterior"]
    args["pc_combination_posterior"] = combinations._replace(
        covariance=combinations.covariance.iloc[::-1, :]
    )
    return args


def test_jet_correction_aligns_original_ids_and_covariance_labels(jets):
    result = API.pc_jet_correction_comparison(**jets)
    assert result.index.equals(jets["groups_index"])
    assert_allclose(result["scale"], [320, 20, 160, 80])
    assert_allclose(result["correction"], -4 / np.array([320, 20, 160, 80]))
    # Var[H(y)] = (1-y)^2 Var[H0] + y^2 Var[H1] + 2y(1-y) Cov[H0,H1].
    y = np.array([1, 0.5, 0.5, 0])
    physical_variance = (1 - y) ** 2 + 8 * y**2 + 1.6 * y * (1 - y)
    assert_allclose(result["correction_std"] ** 2, physical_variance / result["scale"] ** 2)
    assert_allclose(result["experimental_relative_std"], [0.005, 0.2, 0.04, np.nan])
    assert result["exceeds_experiment"].iloc[:3].tolist() == [True, False, False]
    assert pd.isna(result["exceeds_experiment"].iloc[3])
    assert jets["pc_datasets"][0].cuts is not None


def test_factorised_3d_dijets_use_product_mean_and_combination_band():
    args = _inputs(
        "DIJET_3D",
        {
            "ystar": [0.25, 0.75, 0.25, 0.75],
            "yb": [0.25, 0.25, 0.75, 0.75],
            "m_jj": [100, 200, 300, 400],
        },
        {
            "H2j_ystar": {"nodes": [0, 1], "yshift": [2, 2]},
            "H2j_yb": {"nodes": [0, 1], "yshift": [3, 3]},
        },
        np.ones(4),
        0.2 * np.eye(4) + 0.05 * np.ones((4, 4)),
        [10, 10, 10, 10],
        [3, 0, 2, 1],
    )
    result = API.pc_jet_correction_comparison(**args)
    # Interpolated means are 4 and 6; every cross-covariance is 2.4.
    assert_allclose(result["correction"], 26.4 / result["scale"])
    # Prior responses are 6 w_i(ystar) w_j(yb) / m, with sum_i w_i = 1.
    a, b = result["rapidity"], result["yb"]
    sum_squared_weights = (a**2 + (1 - a) ** 2) * (b**2 + (1 - b) ** 2)
    assert_allclose(
        result["correction_std"] ** 2,
        36 * (0.2 * sum_squared_weights + 0.05) / result["scale"] ** 2,
    )
    summary = comparison.pc_jet_correction_exceedance(result)
    assert len(summary) == 4 and summary["ndata"].sum() == 4
    figures = list(comparison.plot_pc_jet_correction_comparison(result))
    try:
        assert len(figures) == 1 and len(figures[0].axes) == 2
        figures[0].canvas.draw()
        assert [ax.get_title() for ax in figures[0].axes] == ["$y_b=0.25$", "$y_b=0.75$"]
    finally:
        for fig in figures:
            plt.close(fig)


@pytest.mark.parametrize(
    "experiment,rapidity,parameter", [("CMS", "ymax", "H2j_ymax"), ("ATLAS", "ystar", "H2j_ystar")]
)
def test_one_dimensional_dijets(experiment, rapidity, parameter):
    args = _inputs(
        "DIJET",
        {rapidity: [0.25, 0.75], "m_jj": [100, 200]},
        {parameter: {"nodes": [0, 1], "yshift": [2, 2]}},
        [1, 1],
        np.eye(2),
        [10, 10],
        [1, 0],
    )
    args["pc_datasets"][0].commondata.metadata.experiment = experiment
    result = comparison.pc_jet_correction_comparison(**args)
    assert result["rapidity_variable"].tolist() == [rapidity] * 2
    assert_allclose(result["correction"], [0.01, 0.02])


def test_exceedance_uses_absolute_correction_and_highest_observed_scale(jets):
    points = comparison.pc_jet_correction_comparison(**jets)
    # Nonmonotonic exceedance; a zero-data point remains undefined.
    points["rapidity"] = 0.5
    points["exceeds_experiment"] = pd.array([True, True, False, pd.NA], dtype="boolean")
    summary = API.pc_jet_correction_exceedance(pc_jet_correction_comparison=points)
    assert_allclose(
        summary[["ndata", "n_compared", "n_exceeding", "max_scale_exceeding"]], [[4, 3, 2, 320]]
    )
    points["exceeds_experiment"] = pd.array([False] * 4, dtype="boolean")
    assert np.isnan(comparison.pc_jet_correction_exceedance(points)["max_scale_exceeding"].iloc[0])


def test_multiple_fits_resolve_and_render_their_own_tables(jets):
    points = comparison.pc_jet_correction_comparison(**jets)
    called = []

    def resolve(fit):
        called.append(fit)
        return [points if fit == "first" else points.iloc[:1]]

    fits = [SimpleNamespace(name=name) for name in ("first", "second")]
    result = comparison.pc_fits_jet_correction_comparison(
        fits=fits, pc_merge_api=SimpleNamespace(pc_fit_jet_correction_comparison=resolve)
    )
    assert called == ["first", "second"]
    assert [len(result.xs(fit.name, level="fit")) for fit in fits] == [4, 1]
    summary = API.pc_fits_jet_correction_exceedance(pc_fits_jet_correction_comparison=result)
    assert summary.groupby(level="fit")["ndata"].sum().tolist() == [4, 1]
    figures = list(comparison.plot_pc_fits_jet_correction_comparison(result))
    try:
        assert [fig.get_suptitle() for fig in figures] == ["first\ntest", "second\ntest"]
        for fig in figures:
            fig.canvas.draw()
    finally:
        for fig in figures:
            plt.close(fig)
    with pytest.raises(CheckError, match="distinct"):
        comparison.pc_fits_jet_correction_comparison([fits[0], fits[0]], None)


def test_no_jet_inputs_return_empty_tables_and_figures(jets):
    jets["pc_datasets"] = []
    empty = API.pc_jet_correction_comparison(**jets)
    assert empty.empty
    assert comparison.pc_jet_correction_exceedance(empty).empty
    assert not list(comparison.plot_pc_jet_correction_comparison(empty))
    multiple = pd.concat({"no_jets": empty}, names=["fit"])
    assert comparison.pc_fits_jet_correction_exceedance(multiple).empty
    assert not list(comparison.plot_pc_fits_jet_correction_comparison(multiple))
