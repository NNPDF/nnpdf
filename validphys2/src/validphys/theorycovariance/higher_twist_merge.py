"""Merge dataset-owned PC determinations, independently or conditional on DIS."""

from collections import namedtuple
import logging

from matplotlib.colors import Normalize, SymLogNorm
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve

from reportengine import collect
from reportengine.checks import check, make_argcheck
from reportengine.figure import figure
from reportengine.table import table
from validphys import plotutils
from validphys.theorycovariance import higher_twist_posterior as pc
from validphys.theorycovariance.output import matrix_plot_labels

PCMergeFit = namedtuple("PCMergeFit", "name parameters basis posterior shifts datasets pdf")
PCMergedPosterior = namedtuple("PCMergedPosterior", "mean covariance owners")
_DIS_PARAMETERS = ("f2p", "f2d", "dis_cc")
log = logging.getLogger(__name__)


def pc_merge_fit(
    fit,
    pc_parameters,
    pc_parameter_basis,
    pc_combination_posterior,
    pc_shifts,
    pc_datasets,
    pc_covmat_pdf,
):
    """Compact per-fit inputs; retain no full experimental or theory matrices."""
    return PCMergeFit(
        fit.name,
        pc_parameters,
        pc_parameter_basis,
        pc_combination_posterior,
        pc_shifts,
        pc_datasets,
        pc_covmat_pdf,
    )


pc_merge_fit_context = collect("pc_merge_fit", ("pc_fit_context",))


def pc_merge_inputs(fits, pc_merge_api):
    """Resolve one fit at a time, releasing its heavy reportengine graph afterwards.

    Collecting all fits in the outer graph would retain every fit's full
    experimental covariance and prediction intermediates simultaneously.
    """
    check(len(fits) > 0, "Provide at least one fit to merge.")
    check(len({fit.name for fit in fits}) == len(fits), "Merge fits must be distinct.")
    inputs = []
    for fit in fits:
        log.info("Extracting PC merge inputs from %s", fit.name)
        inputs.append(pc_merge_api.pc_merge_fit_context(fit=fit.name)[0])
    return inputs


def _reference(inputs, name):
    check(len(inputs) > 0, "Provide at least one fit to merge.")
    if name is None:
        return inputs[0]
    matches = [item for item in inputs if item.name == name]
    check(len(matches) == 1, f"Unknown or repeated reference fit: {name}.")
    return matches[0]


def _dis_basis(item, reference):
    """Align physical DIS nodes, allowing different widths but identical fixed nodes."""
    names = [name for name in _DIS_PARAMETERS if name in reference.parameters]
    check(names, "The conditional merge needs DIS PC parameters in the reference fit.")
    check(
        set(names) == set(_DIS_PARAMETERS).intersection(item.parameters),
        f"DIS parameter types differ between {item.name} and {reference.name}.",
    )
    for name in names:
        par, ref = item.parameters[name], reference.parameters[name]
        check(
            np.array_equal(par["nodes"], ref["nodes"]),
            f"DIS node positions differ for {name} in {item.name}.",
        )
        check(
            np.array_equal(np.asarray(par["yshift"]) == 0, np.asarray(ref["yshift"]) == 0),
            f"Fixed DIS nodes differ for {name} in {item.name}.",
        )
    labels = reference.basis.index[reference.basis.index.get_level_values("HT").isin(names)]
    return item.basis.loc[labels]


def _dis_posterior(item, reference):
    basis = _dis_basis(item, reference)
    return basis, basis @ item.posterior.mean, basis @ item.posterior.covariance @ basis.T


def _conditional_transport(item, reference):
    """Keep p(z_k | theta_D), replacing its DIS marginal by the reference one.

    z_k are the original shift-combination coefficients. This also retains
    the factorised-dijet modes lost by projection onto factor-node parameters.
    Remove exactly fixed DIS nodes before solving the conditional regression.
    """
    basis, mean_d, cov_d = _dis_posterior(item, reference)
    _, ref_mean, ref_cov = _dis_posterior(reference, reference)
    active = np.any(basis.to_numpy() != 0, axis=1)
    basis = basis.to_numpy()[active]
    dd = cov_d.to_numpy()[np.ix_(active, active)]
    ref_dd = ref_cov.to_numpy()[np.ix_(active, active)]
    covariance = item.posterior.covariance.to_numpy()
    cross = covariance @ basis.T
    if active.any():
        scale = np.sqrt(np.diag(dd))
        check(np.all(scale > 0), f"Non-positive DIS posterior variance in {item.name}.")
        try:
            factor = cho_factor(dd / np.outer(scale, scale))
        except np.linalg.LinAlgError as error:
            check(
                False,
                f"DIS covariance in {item.name} is not positive definite on free nodes: {error}",
            )
        gain = cho_solve(factor, (cross / scale).T).T / scale
    else:
        gain = np.zeros((len(covariance), 0))
    residual = covariance - gain @ cross.T
    residual = (residual + residual.T) / 2
    transported = residual + gain @ ref_dd @ gain.T
    mean = item.posterior.mean.to_numpy() + gain @ (ref_mean - mean_d).to_numpy()[active]
    posterior = pc.PCPosterior(
        pd.Series(mean, index=item.posterior.mean.index),
        pd.DataFrame(
            transported, index=item.posterior.mean.index, columns=item.posterior.mean.index
        ),
        item.posterior.prior_covariance,
    )
    return posterior, residual, gain, ref_dd


def _owned_indices(inputs, reference, overrides, selected):
    """First occurrence owns shared datasets, with the reference considered first."""
    ordered = [reference] + [item for item in inputs if item.name != reference.name]
    available = {item.name: set(item.shifts.index.get_level_values("dataset")) for item in ordered}
    owners = {}
    for item in ordered:
        for dataset in item.shifts.index.get_level_values("dataset"):
            owners.setdefault(dataset, item.name)
    check(overrides is None or isinstance(overrides, dict), "pc_dataset_owners must be a mapping.")
    for dataset, owner in (overrides or {}).items():
        check(
            owner in available and dataset in available[owner],
            f"Dataset {dataset} is not available in its requested owner {owner}.",
        )
        owners[dataset] = owner
    if selected is not None:
        check(
            isinstance(selected, list) and len(selected) == len(set(selected)),
            "pc_merge_datasets must be a list of distinct dataset names.",
        )
        check(set(selected) <= owners.keys(), "pc_merge_datasets contains unavailable datasets.")
        owners = {name: owners[name] for name in selected}
    result = []
    for item in ordered:
        index = item.shifts.index
        mask = [owners.get(name) == item.name for name in index.get_level_values("dataset")]
        index = index[mask]
        check(index.droplevel(0).is_unique, f"Duplicate data-point labels in {item.name}.")
        if len(index):
            result.append((item, index))
    check(result, "No owned data points selected for the PC merge.")
    return result


def _check_merge_method(pc_merge_method):
    check(
        pc_merge_method in ("independent", "conditional"),
        "pc_merge_method must be 'independent' or 'conditional'.",
    )


@make_argcheck(_check_merge_method)
def pc_merged_posterior(
    pc_merge_inputs,
    pc_merge_method="independent",
    pc_merge_reference=None,
    pc_dataset_owners=None,
    pc_merge_datasets=None,
):
    """Owned mean shifts and Cpost, before covariance inflation and the mean outer product.

    Defaults: the first fit anchors DIS and owns every shared dataset. Optional
    ownership overrides map dataset names to fit names. pc_merge_datasets
    restricts the output; otherwise the union of fitted datasets is included.
    Only the owner's fitted points are retained, with their original point IDs.
    """
    _check_merge_method(pc_merge_method)
    reference = _reference(pc_merge_inputs, pc_merge_reference)
    owned = _owned_indices(pc_merge_inputs, reference, pc_dataset_owners, pc_merge_datasets)
    means, owners, blocks, common = [], [], [], []
    for item, index in owned:
        beta = item.shifts.loc[index].to_numpy()
        if pc_merge_method == "conditional":
            check(
                item.name == reference.name
                or not any(str(group).startswith("DIS") for group in index.get_level_values(0)),
                "In a conditional merge all DIS datasets must be owned by the reference fit.",
            )
            posterior, residual, gain, ref_dd = _conditional_transport(item, reference)
            common.append(beta @ gain)
        else:
            posterior = item.posterior
            residual = posterior.covariance.to_numpy()
        nodes = pc.pc_posterior(posterior, item.basis)
        means.append(pc._prediction_mean(item.datasets, item.pdf, item.parameters, nodes, index))
        owners.append(pd.Series(item.name, index=index, name="owner"))
        blocks.append(beta @ residual @ beta.T)
    mean, owner = pd.concat(means), pd.concat(owners)
    if common:
        shared = np.vstack(common)
        covariance = shared @ ref_dd @ shared.T
    else:
        covariance = np.zeros((len(mean), len(mean)))
    start = 0
    for block in blocks:
        end = start + len(block)
        covariance[start:end, start:end] += block
        start = end
    covariance = (covariance + covariance.T) / 2
    return PCMergedPosterior(
        mean, pd.DataFrame(covariance, index=mean.index, columns=mean.index), owner
    )


@table
def pc_merged_prediction_covmat(pc_merged_posterior):
    """Export the merged Cpost, without the factor or the mean outer product."""
    return pc_merged_posterior.covariance


@table
def pc_merged_mean_shift(pc_merged_posterior):
    """Owned mean shifts and source fit, aligned with the merged matrix."""
    return pd.concat([pc_merged_posterior.mean, pc_merged_posterior.owners], axis=1)


@table
@pc.check_pc_covariance_factor
def pc_merged_covmat(pc_merged_posterior, pc_covariance_factor=2.0):
    """Export factor * Cpost_merged + dT dT.T; the factor scales only Cpost."""
    return pc.pc_shift_covmat(
        pc_merged_posterior.covariance, pc_merged_posterior.mean, pc_covariance_factor
    )


def _plot_merged_matrix(matrix, title, norm, colorbar_label):
    """Show nonzero-variance rows in export order, with dataset labels and boundaries."""
    active = np.diag(matrix) > 0
    # An identically zero covariance is valid (e.g. f=0 and zero mean shifts).
    if active.any():
        matrix = matrix.loc[active, active]
    fig, ax = plotutils.subplots(figsize=(15, 13))
    image = ax.matshow(matrix, cmap="Spectral_r", norm=norm)
    fig.colorbar(image, ax=ax, label=colorbar_label)
    ticks, names, boundaries = matrix_plot_labels(matrix.droplevel(0))
    ax.set_xticks(ticks, names, rotation=90, fontsize=7)
    ax.set_yticks(ticks, names, fontsize=7)
    ax.xaxis.tick_bottom()
    for boundary in boundaries[1:-1]:
        ax.axvline(boundary - 0.5, color="grey", linestyle="--", linewidth=0.5)
        ax.axhline(boundary - 0.5, color="grey", linestyle="--", linewidth=0.5)
    ax.set(title=title, xlabel="Dataset", ylabel="Dataset")
    return fig


@figure
def plot_pc_merged_covmat(pc_merged_covmat):
    """Final f*Cpost + dT dT.T in prediction-unit products, on a symmetric log scale."""
    limit = np.abs(pc_merged_covmat.to_numpy()).max() or 1.0
    return _plot_merged_matrix(
        pc_merged_covmat,
        "Merged PC covariance",
        SymLogNorm(linthresh=limit * 1e-6, vmin=-limit, vmax=limit),
        "Covariance [product of prediction units]",
    )


@figure
def plot_pc_merged_corrmat(pc_merged_covmat):
    """Correlations of the final merged error model, including the mean outer product."""
    return _plot_merged_matrix(
        pc._correlation(pc_merged_covmat),
        "Merged PC correlations",
        Normalize(vmin=-1, vmax=1),
        "Correlation",
    )


@table
def pc_merge_dis_comparison(pc_merge_inputs, pc_merge_reference=None):
    """Compare the input DIS node means and widths before conditional transport.

    Differences in units of combined marginal errors are descriptive, not
    significance tests: all input fits share DIS measurements.
    """
    reference = _reference(pc_merge_inputs, pc_merge_reference)
    _, ref_mean, ref_cov = _dis_posterior(reference, reference)
    ref_std = np.sqrt(np.maximum(0, np.diag(ref_cov)))
    frames = {}
    for item in pc_merge_inputs:
        _, mean, covariance = _dis_posterior(item, reference)
        std = np.sqrt(np.maximum(0, np.diag(covariance)))
        scale = np.sqrt(std**2 + ref_std**2)
        difference = (mean - ref_mean).to_numpy()
        frames[item.name] = pd.DataFrame(
            {
                "mean": mean,
                "std": std,
                "reference_mean": ref_mean,
                "reference_std": ref_std,
                "difference_over_combined_std": np.divide(
                    difference, scale, out=np.zeros_like(difference), where=scale != 0
                ),
            }
        )
    return pd.concat(frames, names=["fit"])


@table
def pc_merge_dis_covmats(pc_merge_inputs, pc_merge_reference=None):
    """Input DIS covariance matrices, stacked by fit for compatibility checks."""
    reference = _reference(pc_merge_inputs, pc_merge_reference)
    return pd.concat(
        {item.name: _dis_posterior(item, reference)[2] for item in pc_merge_inputs}, names=["fit"]
    )
