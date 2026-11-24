"""Power-correction posteriors, evaluated in the ``pc_fit_context`` namespace.

The coefficient-space update follows Eqs. (3.37–3.38) of arXiv:2105.05114.
The prior is expressed through independent shift combinations. PDF fluctuations
are projected into coefficient space without constructing their data covariance.
"""

from collections import namedtuple
from copy import copy
from numbers import Real
import re

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve, solve

from reportengine.checks import check, make_argcheck
from reportengine.figure import figure, figuregen
from reportengine.table import table
from validphys import plotutils
from validphys.tableloader import parse_exp_mat
from validphys.theorycovariance.higher_twist_functions import (
    compute_deltas_pc,
    construct_pars_combs,
    linear_bin_function,
)
from validphys.theorycovariance.output import matrix_plot_labels, plot_corrmat_heatmap

PCPosterior = namedtuple("PCPosterior", "mean covariance prior_covariance")
_DIJET_FACTORS = ("H2j_ystar", "H2j_yb")


def pc_parameter_basis(pc_parameters):
    """Map the prior's independent shift combinations to named node parameters.

    This is diagonal for 1D corrections. For factorised dijets, each product
    combination contributes to both participating nodes.
    The resulting prior can be singular; no inversion in parameter space is used.
    """
    pair = _DIJET_FACTORS if set(_DIJET_FACTORS) <= pc_parameters.keys() else None
    combinations = construct_pars_combs(pc_parameters, pair)
    index = pd.MultiIndex.from_tuples(
        [
            (key, f"{key}({i})")
            for key, par in pc_parameters.items()
            for i in range(len(par["nodes"]))
        ],
        names=["HT", "nodes"],
    )
    return pd.DataFrame(
        np.array([np.concatenate(list(c["comb"].values())) for c in combinations]).T,
        index=index,
        columns=[c["label"] for c in combinations],
    )


def pc_datasets(groups_data_by_process, pc_included_procs, pc_excluded_datasets):
    """Datasets to which the fit applied the power-correction prior."""
    return [
        ds
        for group in groups_data_by_process
        if group.name in pc_included_procs
        for ds in group.datasets
        if ds.name not in pc_excluded_datasets
    ]


def pc_shifts(pc_datasets, pc_covmat_pdf, pc_parameters, pc_parameter_basis, groups_index):
    """Prior shifts (beta), with zero rows for datasets without corrections."""
    beta = pd.DataFrame(0.0, index=groups_index, columns=pc_parameter_basis.columns)
    for ds in pc_datasets:
        shifts = compute_deltas_pc(ds, pc_covmat_pdf, pc_parameters)
        rows = beta.index.get_level_values("dataset") == ds.name
        beta.loc[rows, :] = np.column_stack([shifts[key] for key in beta.columns])
    return beta


def pc_central_pseudodata(
    read_pdf_pseudodata, fit_diagonal_basis_rotation, data_index, groups_index
):
    """Mean postfit pseudodata in prediction order, undoing the fitting rotation."""
    check(len(read_pdf_pseudodata) > 0, "No saved PDF pseudodata replicas found.")
    rotation = fit_diagonal_basis_rotation
    if rotation is None:
        target = groups_index.droplevel(0)
        samples = [
            p.pseudodata.droplevel(0).loc[target].to_numpy()[:, 0] for p in read_pdf_pseudodata
        ]
        central = np.mean(samples, axis=0)
    else:
        modes = sorted(read_pdf_pseudodata[0].pseudodata.index, key=lambda x: int(x.split()[-1]))
        check(len(modes) == len(data_index), "Pseudodata and fitting rotation sizes differ.")
        central = np.mean(
            [p.pseudodata.loc[modes].to_numpy()[:, 0] for p in read_pdf_pseudodata], axis=0
        )
        central = pd.Series(solve(rotation, central), index=data_index.droplevel(0))
        central = central.loc[groups_index.droplevel(0)].to_numpy()
    check(np.isfinite(central).all(), "Non-finite central pseudodata.")
    return pd.Series(central, index=groups_index)


def _read_covmat(path, index):
    """Align both axes by (dataset, id), independently of the plotting grouping."""
    matrix = parse_exp_mat(path).droplevel(0).droplevel(0, axis=1)
    target = index.droplevel(0)
    check(matrix.index.is_unique and matrix.columns.is_unique, f"Duplicate data labels in {path}.")
    matrix = matrix.loc[target, target].to_numpy()
    check(np.isfinite(matrix).all(), f"Non-finite theory covariance in {path}.")
    return matrix


def _condition(beta, basis, total_covariance, predictions, central_data):
    """Gaussian update with population (1/Nrep) PDF covariance."""
    prior = basis @ basis.T
    cross = basis @ beta.T
    gain = cho_solve(cho_factor(total_covariance), cross.T).T
    central = predictions.mean(axis=1)
    projected = gain @ (predictions - central[:, None])
    mean = gain @ (central_data - central)
    covariance = prior - gain @ cross.T + projected @ projected.T / predictions.shape[1]
    return mean, (covariance + covariance.T) / 2, prior


def pc_combination_posterior(
    fit,
    pc_shifts,
    experiments_covmat_no_table,
    group_result_table_no_table,
    pc_central_pseudodata,
    pdfreplicas,
):
    """Posterior of the unit-variance prior combinations, including PDF uncertainty.

    Read the stored prescription matrices one at a time. The experimental
    covariance uses data-normalised multiplicative errors.
    Build it by experiment to retain correlations across process groups.
    All non-PC theory contributions remain in the conditioning covariance.
    """
    index = pc_shifts.index
    beta = pc_shifts.to_numpy()
    experimental = experiments_covmat_no_table.droplevel(0).droplevel(0, axis=1)
    target = index.droplevel(0)
    total = experimental.loc[target, target].to_numpy(copy=True)
    settings = fit.as_input()["theorycovmatconfig"]
    for i, prescription in enumerate(settings["point_prescriptions"]):
        path = (
            fit.path
            / "tables"
            / (
                f"datacuts_theory_theorycovmatconfig_point_prescriptions{i}"
                "_theory_covmat_custom_per_prescription.csv"
            )
        )
        covariance = _read_covmat(path, index)
        if prescription == "power corrections":
            # Check in blocks to avoid retaining a second Ndata x Ndata matrix.
            for start in range(0, len(index), 256):
                rows = slice(start, start + 256)
                check(
                    np.allclose(beta[rows] @ beta.T, covariance[rows], rtol=1e-4, atol=1e-8),
                    "Reconstructed PC shifts disagree with the fit's saved covariance. "
                    "Check the prior PDF, dataset exclusions and PC parametrisation.",
                )
        total += covariance
    if settings.get("user_covmat_path") is not None:
        total += _read_covmat(
            fit.path / "tables" / "datacuts_theory_theorycovmatconfig_user_covmat.csv", index
        )
    predictions = group_result_table_no_table.loc[index].filter(regex=r"^rep_\d+$").to_numpy()
    check(
        predictions.shape[1] == len(pdfreplicas) and len(pdfreplicas) > 0,
        "Predictions and saved pseudodata must use the same postfit replicas.",
    )
    mean, covariance, prior = _condition(
        beta, np.eye(beta.shape[1]), total, predictions, pc_central_pseudodata.loc[index].to_numpy()
    )
    parameters = pc_shifts.columns
    return PCPosterior(
        pd.Series(mean, index=parameters, name="mean"),
        pd.DataFrame(covariance, index=parameters, columns=parameters),
        pd.DataFrame(prior, index=parameters, columns=parameters),
    )


def pc_posterior(pc_combination_posterior, pc_parameter_basis):
    """Project the common combination posterior onto the named node parameters."""
    basis = pc_parameter_basis
    return PCPosterior(
        basis @ pc_combination_posterior.mean,
        basis @ pc_combination_posterior.covariance @ basis.T,
        basis @ basis.T,
    )


@table
def pc_posterior_table(pc_posterior, pc_parameters):
    """Node locations, posterior means and marginal prior/posterior uncertainties."""
    return pd.DataFrame(
        {
            "node": np.concatenate([par["nodes"] for par in pc_parameters.values()]),
            "mean": pc_posterior.mean,
            "std": np.sqrt(np.maximum(0, np.diag(pc_posterior.covariance))),
            "prior_std": np.sqrt(np.diag(pc_posterior.prior_covariance)),
        }
    )


@table
def pc_posterior_covmat(pc_posterior):
    """Full parameter covariance, including correlations between PC types."""
    return pc_posterior.covariance


@table
def pc_posterior_corrmat(pc_posterior):
    """Parameter correlations; fixed nodes have zero rows and columns."""
    return _correlation(pc_posterior.covariance)


def _correlation(covariance):
    """Normalise a labelled covariance, assigning zero to fixed-point rows."""
    std = np.sqrt(np.maximum(0, np.diag(covariance)))
    norm = np.outer(std, std)
    correlation = np.divide(covariance, norm, out=np.zeros_like(norm), where=norm != 0)
    return pd.DataFrame(correlation, index=covariance.index, columns=covariance.columns)


def pc_prediction_covmat(pc_combination_posterior, pc_shifts):
    """Posterior covariance of PC shifts at fitted points: B P_combination B.T.

    This is the centred covariance in the fitted shift-combination model,
    including correlations across datasets. It also supports singular node
    bases (factorised dijets), without inverting the node-to-combination map.
    """
    return pc_shifts @ pc_combination_posterior.covariance @ pc_shifts.T


def _prediction_mean(datasets, pdf, parameters, posterior, index):
    """Evaluate the correction mean from the node posterior.

    Products include their cross-covariance. Ratio corrections are evaluated
    at the node means, rather than taking the expectation of a Gaussian ratio.
    """
    shifted = {
        name: {**par, "yshift": posterior.mean.xs(name, level="HT").to_numpy()}
        for name, par in parameters.items()
    }
    factorised = set(_DIJET_FACTORS) <= parameters.keys()
    unit = {name: {**par, "yshift": np.ones(len(par["nodes"]))} for name, par in parameters.items()}
    result = pd.Series(0.0, index=index, name="mean_shift")
    for ds in datasets:
        rows = index.get_level_values("dataset") == ds.name
        if not rows.any():
            continue
        values = compute_deltas_pc(ds, pdf, shifted, _all_nodes)["all"]
        if factorised:
            responses = compute_deltas_pc(ds, pdf, unit)
            a, b = _DIJET_FACTORS
            for i in range(len(parameters[a]["nodes"])):
                for j in range(len(parameters[b]["nodes"])):
                    values = (
                        values
                        + responses[f"{a}({i})_{b}({j})"]
                        * posterior.covariance.loc[(a, f"{a}({i})"), (b, f"{b}({j})")]
                    )
        result.loc[rows] = values
    check(np.isfinite(result).all(), "Non-finite PC mean shifts.")
    return result


def pc_prediction_mean(pc_datasets, pc_covmat_pdf, pc_parameters, pc_posterior, groups_index):
    """Per-fit mean PC shift, including the factorised-product covariance term."""
    return _prediction_mean(pc_datasets, pc_covmat_pdf, pc_parameters, pc_posterior, groups_index)


@table
def pc_prediction_mean_table(pc_prediction_mean):
    """Export the mean PC shift in prediction units."""
    return pc_prediction_mean.to_frame()


@make_argcheck
def check_pc_covariance_factor(pc_covariance_factor):
    check(
        isinstance(pc_covariance_factor, Real)
        and not isinstance(pc_covariance_factor, bool)
        and np.isfinite(pc_covariance_factor)
        and pc_covariance_factor >= 0,
        "pc_covariance_factor must be a finite, non-negative number.",
    )


@table
@check_pc_covariance_factor
def pc_shift_covmat(pc_prediction_covmat, pc_prediction_mean, pc_covariance_factor=2.0):
    """Inflated second moment: factor * Cpost + mean_shift mean_shift.T."""
    mean = pc_prediction_mean.loc[pc_prediction_covmat.index].to_numpy()
    return pc_covariance_factor * pc_prediction_covmat + np.outer(mean, mean)


@table
def pc_prediction_covmat_table(pc_prediction_covmat):
    """Export the full prediction-space PC posterior covariance."""
    return pc_prediction_covmat


def pc_prediction_corrmat(pc_prediction_covmat):
    """Prediction-space PC correlations; unaffected points have zero rows."""
    return _correlation(pc_prediction_covmat)


@table
def pc_prediction_corrmat_table(pc_prediction_corrmat):
    """Export the full prediction-space PC posterior correlation matrix."""
    return pc_prediction_corrmat


@table
def pc_prediction_uncertainties(pc_combination_posterior, pc_shifts, group_result_table_no_table):
    """Prior/posterior PC standard deviations, also relative to the mean prediction.

    The diagonal is evaluated directly, without building a data-space matrix.
    Relative uncertainties are undefined (NaN) where the mean prediction is zero.
    """
    beta = pc_shifts.to_numpy()
    posterior = pc_combination_posterior.covariance.loc[pc_shifts.columns, pc_shifts.columns]
    prior_std = np.sqrt(np.einsum("ij,ij->i", beta, beta))
    std = np.sqrt(np.maximum(0, np.einsum("ij,jk,ik->i", beta, posterior, beta)))
    central = (
        group_result_table_no_table.loc[pc_shifts.index].filter(regex=r"^rep_\d+$").mean(axis=1)
    )
    scale = np.abs(central.to_numpy())
    scale[scale == 0] = np.nan
    return pd.DataFrame(
        {
            "theory_mean": central,
            "prior_std": prior_std,
            "posterior_std": std,
            "prior_relative_std": prior_std / scale,
            "posterior_relative_std": std / scale,
        },
        index=pc_shifts.index,
    )


@table
def pc_prediction_experimental_uncertainties(
    pc_prediction_uncertainties, experiments_covmat_no_table, group_result_table_no_table
):
    """Compare sqrt(diag(Cpost)) with experimental errors at the same fitted points.

    Both relative errors use the absolute data central value.
    Cpost is uninflated and excludes the mean outer product.
    Zero denominators give NaN. Align by (dataset, id), since the experimental
    covariance is grouped by experiment rather than by process.
    """
    points = pc_prediction_uncertainties
    index = points.index
    experimental = pd.Series(
        np.sqrt(np.diag(experiments_covmat_no_table)),
        index=experiments_covmat_no_table.index.droplevel(0),
    ).loc[index.droplevel(0)]
    experimental = pd.Series(experimental.to_numpy(), index=index)
    central = group_result_table_no_table.loc[index, "data_central"]
    scale = central.abs().replace(0, np.nan)
    return pd.DataFrame(
        {
            "data_central": central,
            "prior_std": points["prior_std"],
            "posterior_std": points["posterior_std"],
            "experimental_std": experimental,
            "posterior_relative_std": points["posterior_std"] / scale,
            "experimental_relative_std": experimental / scale,
            "posterior_over_experimental": points["posterior_std"]
            / experimental.replace(0, np.nan),
        },
        index=index,
    )


@figure
def plot_pc_prediction_corrmat(pc_prediction_corrmat, pc_shifts):
    """PC correlations between fitted points with a nonzero prior correction."""
    active = (pc_shifts != 0).any(axis=1)
    check(active.any(), "No fitted points have a nonzero PC prior.")
    return plot_corrmat_heatmap(
        pc_prediction_corrmat.loc[active, active], "Posterior PC correlations in prediction space"
    )


def _plot_prediction_uncertainties(uncertainties, columns, ylabel):
    """Pointwise comparisons by process, with the same PC selection and dataset labels."""
    active = uncertainties.query("prior_std > 0")
    for process, points in active.groupby(level=0, sort=False):
        fig, ax = plotutils.subplots(figsize=(14, 5))
        for column, label in columns:
            ax.plot(100 * points[column].to_numpy(), ".", markersize=3, label=label)
        ticks, labels, boundaries = matrix_plot_labels(points.droplevel(0))
        ax.set_xticks(ticks, labels, rotation=90, fontsize=8)
        for boundary in boundaries[1:-1]:
            ax.axvline(boundary - 0.5, color="grey", linestyle="--", linewidth=0.5)
        ax.set(title=process, ylabel=ylabel, ylim=(0, None))
        ax.margins(x=0)
        ax.legend()
        yield fig


@figuregen
def plot_pc_prediction_uncertainties(pc_prediction_uncertainties):
    """Prior and posterior relative PC uncertainties, grouped by process."""
    yield from _plot_prediction_uncertainties(
        pc_prediction_uncertainties,
        (("prior_relative_std", "Prior PC"), ("posterior_relative_std", "Posterior PC")),
        r"PC uncertainty / $|\overline{T}|$ [%]",
    )


@figuregen
def plot_pc_prediction_experimental_uncertainties(pc_prediction_experimental_uncertainties):
    """Posterior PC and experimental errors relative to the absolute data value, by process."""
    yield from _plot_prediction_uncertainties(
        pc_prediction_experimental_uncertainties,
        (("posterior_relative_std", "Posterior PC"), ("experimental_relative_std", "Experiment")),
        r"Uncertainty / $|D|$ [%]",
    )


def _curve_moments(weights, mean, covariance, average=None):
    central = weights @ mean
    variance = np.einsum("ij,jk,ik->i", weights, covariance, weights)
    if average is not None:
        other = average @ mean
        other_var = average @ covariance @ average
        cross = weights @ covariance @ average
        variance = (
            central**2 * other_var
            + other**2 * variance
            + 2 * central * other * cross
            + variance * other_var
            + cross**2
        )
        central = central * other + cross
    return central, np.sqrt(np.maximum(0, variance))


@figuregen
def plot_pc_posterior(pc_posterior, pc_parameters):
    """Prior and posterior one-sigma bands, propagating the full node covariance.

    Factorised dijets are shown as rapidity-averaged corrections in GeV.
    Gaussian product moments are evaluated analytically, without sampling.
    """
    index = pc_posterior.mean.index.get_level_values("HT")
    factorised = set(_DIJET_FACTORS) <= pc_parameters.keys()
    labels = {"Hj": "$y$", "H2j_ystar": "$y^*$", "H2j_yb": "$y_b$", "H2j_ymax": r"$|y|_{\max}$"}
    for name, parameters in pc_parameters.items():
        nodes = np.asarray(parameters["nodes"])
        x = np.unique(np.r_[nodes, np.linspace(nodes[0], nodes[-1], 200)])
        weights = np.zeros((len(x), len(index)))
        weights[:, index == name] = np.array(
            [linear_bin_function(x, row, nodes) for row in np.eye(len(nodes))]
        ).T
        average = None
        if factorised and name in _DIJET_FACTORS:
            other = _DIJET_FACTORS[1 - _DIJET_FACTORS.index(name)]
            other_nodes = np.asarray(pc_parameters[other]["nodes"])
            intervals = np.diff(other_nodes)
            average = np.zeros(len(index))
            average[index == other] = (np.r_[intervals, 0] + np.r_[0, intervals]) / (
                2 * (other_nodes[-1] - other_nodes[0])
            )
        fig, ax = plotutils.subplots()
        for label, mean, cov in (
            ("Prior", np.zeros(len(index)), pc_posterior.prior_covariance.to_numpy()),
            ("Posterior", pc_posterior.mean.to_numpy(), pc_posterior.covariance.to_numpy()),
        ):
            central, std = _curve_moments(weights, mean, cov, average)
            (line,) = ax.plot(x, central, label=label)
            ax.fill_between(x, central - std, central + std, color=line.get_color(), alpha=0.25)
        dis = name in {"f2p", "f2d", "dis_cc"}
        ax.set(
            xlabel="$x$" if dis else labels.get(name, "rapidity"),
            ylabel=r"$H\;[\mathrm{GeV}^2]$" if dis else r"$H\;[\mathrm{GeV}]$",
            title=name,
        )
        ax.legend()
        yield fig


@figure
def plot_pc_posterior_corrmat(pc_posterior_corrmat):
    """Correlation heatmap for all PC/HT node parameters."""
    fig, ax = plotutils.subplots(figsize=(10, 8))
    image = ax.imshow(pc_posterior_corrmat, vmin=-1, vmax=1, cmap="RdBu_r")
    ticks = np.arange(len(pc_posterior_corrmat))
    labels = pc_posterior_corrmat.index.get_level_values("nodes")
    ax.set_xticks(ticks, labels, rotation=90, fontsize=7)
    ax.set_yticks(ticks, labels, fontsize=7)
    fig.colorbar(image, ax=ax, label="Correlation")
    return fig


def _all_nodes(parameters, _pair=None):
    return [{"label": "all", "comb": {key: par["yshift"] for key, par in parameters.items()}}]


def _posterior_cfactor(dataset, parameters, posterior):
    """Full-kinematics c-factors; ratio operands remain separate."""
    dataset = copy(dataset)
    dataset.cuts = None
    shifted = {
        name: {**par, "yshift": posterior.mean.xs(name, level="HT").to_numpy()}
        for name, par in parameters.items()
    }
    result = compute_deltas_pc(dataset, None, shifted, _all_nodes, cfactors=True)["all"]
    if set(_DIJET_FACTORS) <= parameters.keys() and not isinstance(result, tuple):
        unit = {
            name: {**par, "yshift": np.ones(len(par["nodes"]))} for name, par in parameters.items()
        }
        responses = compute_deltas_pc(dataset, None, unit, cfactors=True)
        a, b = _DIJET_FACTORS
        for i in range(len(parameters[a]["nodes"])):
            for j in range(len(parameters[b]["nodes"])):
                covariance = posterior.covariance.loc[(a, f"{a}({i})"), (b, f"{b}({j})")]
                result += covariance * (responses[f"{a}({i})_{b}({j})"] - 1)
    return result if isinstance(result, tuple) else (result,)


@table
def pc_posterior_cfactors(
    fit, pc_datasets, pc_parameters, pc_posterior, output_path, pc_cfactor_label: str = "PC"
):
    """Write central c-factors per FK operand and return a file manifest.

    Files cover all FK bins, including points removed by the fit cuts. The
    uncertainty column is zero: correlated errors are in ``pc_posterior_covmat``.
    """
    from pineappl.fk_table import FkTable

    check(
        re.fullmatch(r"[A-Za-z0-9_-]+", pc_cfactor_label) is not None,
        "pc_cfactor_label must contain only letters, numbers, underscores or hyphens.",
    )
    folder = output_path / "cfactors" / fit.name
    folder.mkdir(parents=True, exist_ok=True)
    records, written = [], {}
    for dataset in pc_datasets:
        factors = _posterior_cfactor(dataset, pc_parameters, pc_posterior)
        check(len(factors) == len(dataset.fkspecs), f"Unsupported FK operands for {dataset.name}.")
        for spec, values in zip(dataset.fkspecs, factors):
            check(not spec.legacy, "PC c-factor export requires PineAPPL FK tables.")
            offset = 0
            for path in spec.fkpath:
                name = path.name.removesuffix(".pineappl.lz4")
                offset += (spec.metadata.shifts or {}).get(name, 0)
                size = len(FkTable.read(path).bin_normalizations())
                check(
                    0 <= offset and offset + size <= len(values),
                    f"FK bins do not match commondata for {name}.",
                )
                cfac = values[offset : offset + size]
                check(np.isfinite(cfac).all(), f"Non-finite c-factors for {name}.")
                filename = folder / f"CF_{pc_cfactor_label}_{name}.dat"
                if filename in written:
                    check(
                        np.allclose(written[filename], cfac), f"Conflicting c-factors for {name}."
                    )
                else:
                    delimiter = "*" * 80
                    header = (
                        f"{delimiter}\nSetName: {dataset.name}\nFromFit: {fit.name}\n"
                        "Description: Posterior power correction\n"
                        f"Uncertainties: pc_posterior_covmat\n{delimiter}"
                    )
                    np.savetxt(
                        filename,
                        np.column_stack((cfac, np.zeros(size))),
                        header=header,
                        comments="",
                        fmt="%.12e",
                    )
                    written[filename] = cfac
                records.append(
                    {
                        "dataset": dataset.name,
                        "fktable": name,
                        "ndata": size,
                        "path": str(filename.relative_to(output_path)),
                    }
                )
                offset += size
    return pd.DataFrame(records)
