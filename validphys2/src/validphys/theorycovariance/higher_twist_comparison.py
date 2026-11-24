"""Compare fitted jet PC c-factors with experimental errors to assess kinematic cuts."""

from copy import copy

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.checks import check
from reportengine.figure import figuregen
from reportengine.table import table
from validphys import plotutils
from validphys.theorycovariance import higher_twist_posterior as pc
from validphys.theorycovariance.higher_twist_functions import DIJET_RAPIDITY_VAR, compute_deltas_pc


@table
def pc_jet_correction_comparison(
    pc_datasets,
    pc_parameters,
    pc_posterior,
    pc_combination_posterior,
    experiments_covmat_no_table,
    groups_index,
):
    """Central C_PC - 1 versus sigma_exp / abs(data), on the selected fit's jet points.

    C-factors are evaluated from the posterior, not read from installed files.
    Their bands propagate the posterior shift-combination covariance in
    dimensionless c-factor units. The exceedance criterion compares the
    absolute central correction with the relative experimental error.
    No inflation or mean outer product enters this comparison.
    """
    experimental = pd.Series(
        np.sqrt(np.diag(experiments_covmat_no_table)),
        index=experiments_covmat_no_table.index.droplevel(0),
    )
    columns = [
        "scale_variable",
        "scale",
        "rapidity_variable",
        "rapidity",
        "yb",
        "correction",
        "correction_std",
        "experimental_relative_std",
        "exceeds_experiment",
    ]
    frames = []
    for dataset in pc_datasets:
        metadata = dataset.commondata.metadata
        process = metadata.process_type.name
        if process not in ("JET", "DIJET", "DIJET_3D"):
            continue
        index = groups_index[groups_index.get_level_values("dataset") == dataset.name]
        if not len(index):
            continue
        ids = index.get_level_values("id").to_numpy(dtype=int)
        kinematics = metadata.load_kinematics()
        scale = "pT" if process == "JET" else "m_jj"
        rapidity = (
            "y"
            if process == "JET"
            else "ystar" if process == "DIJET_3D" else DIJET_RAPIDITY_VAR[metadata.experiment]
        )

        def values(key):
            return kinematics[key].to_numpy().reshape(-1)[ids]

        (correction,) = pc._posterior_cfactor(dataset, pc_parameters, pc_posterior)
        uncut = copy(dataset)
        uncut.cuts = None
        shifts = compute_deltas_pc(uncut, None, pc_parameters, cfactors=True)
        covariance = pc_combination_posterior.covariance
        covariance = covariance.loc[covariance.index, covariance.index]
        response = np.column_stack([shifts[key][ids] - 1 for key in covariance.index])
        variance = np.einsum("ij,jk,ik->i", response, covariance, response)
        central = dataset.commondata.load().central_values.to_numpy(dtype=float)[ids]
        norm = np.abs(central)
        norm[norm == 0] = np.nan
        frame = pd.DataFrame(
            {
                "scale_variable": scale,
                "scale": values(scale),
                "rapidity_variable": rapidity,
                "rapidity": values(rapidity),
                "yb": values("yb") if process == "DIJET_3D" else np.nan,
                "correction": correction[ids] - 1,
                "correction_std": np.sqrt(np.maximum(0, variance)),
                "experimental_relative_std": experimental.loc[index.droplevel(0)].to_numpy() / norm,
            },
            index=index,
        )
        frame["exceeds_experiment"] = (
            (frame["correction"].abs() > frame["experimental_relative_std"])
            .astype("boolean")
            .where(frame["experimental_relative_std"].notna())
        )
        frames.append(frame)
    return pd.concat(frames) if frames else pd.DataFrame(columns=columns, index=groups_index[:0])


@table
def pc_jet_correction_exceedance(pc_jet_correction_comparison):
    """Highest observed pT or m_jj satisfying abs(C_PC - 1) > sigma_exp / abs(data).

    Report each rapidity bin separately (yb, ystar for 3D dijets). NaN means
    no compared point exceeds the error. This does not assume monotonicity or
    extrapolate beyond the owner's fitted points.
    """
    records = []
    for dataset, points in pc_jet_correction_comparison.groupby(level="dataset", sort=False):
        for (yb, rapidity), bins in points.groupby(["yb", "rapidity"], sort=True, dropna=False):
            valid = bins["exceeds_experiment"].notna()
            exceeding = bins["exceeds_experiment"].fillna(False)
            records.append(
                {
                    "dataset": dataset,
                    "yb": yb,
                    "rapidity": rapidity,
                    "rapidity_variable": bins["rapidity_variable"].iloc[0],
                    "scale_variable": bins["scale_variable"].iloc[0],
                    "ndata": len(bins),
                    "n_compared": int(valid.sum()),
                    "n_exceeding": int(exceeding.sum()),
                    "max_scale_exceeding": bins.loc[exceeding, "scale"].max(),
                    "min_scale": bins["scale"].min(),
                    "max_scale": bins["scale"].max(),
                    "max_abs_correction": bins["correction"].abs().max(),
                }
            )
    columns = [
        "dataset",
        "yb",
        "rapidity",
        "rapidity_variable",
        "scale_variable",
        "ndata",
        "n_compared",
        "n_exceeding",
        "max_scale_exceeding",
        "min_scale",
        "max_scale",
        "max_abs_correction",
    ]
    return pd.DataFrame(records, columns=columns).set_index(["dataset", "yb", "rapidity"])


@figuregen
def plot_pc_jet_correction_comparison(pc_jet_correction_comparison):
    """Jet corrections with posterior bands, experimental errors and rapidity-bin separators."""
    rapidity_labels = {"y": "|y|", "ystar": "y^*", "ymax": r"|y|_{\max}"}
    for dataset, points in pc_jet_correction_comparison.groupby(level="dataset", sort=False):
        outer = [(yb, row) for yb, row in points.groupby("yb", sort=True, dropna=False)]
        fig, axes = plotutils.subplots(
            nrows=len(outer), figsize=(15, 4 * len(outer)), squeeze=False
        )
        fig.suptitle(dataset)
        for ax, (yb, row) in zip(axes[:, 0], outer):
            offset, ticks, labels = 0, [], []
            for n, (rapidity, bins) in enumerate(row.groupby("rapidity", sort=True)):
                bins = bins.sort_values("scale")
                x = np.arange(len(bins)) + offset
                mean = 100 * bins["correction"].to_numpy()
                std = 100 * bins["correction_std"].to_numpy()
                (line,) = ax.plot(
                    x, mean, ".", color="C0", label="Mean PC correction" if n == 0 else None
                )
                ax.fill_between(
                    x,
                    mean - std,
                    mean + std,
                    color=line.get_color(),
                    alpha=0.25,
                    label="Posterior 1-sigma band" if n == 0 else None,
                )
                ax.plot(
                    x,
                    100 * bins["experimental_relative_std"],
                    ".",
                    color="C1",
                    label="Experimental relative error" if n == 0 else None,
                )
                if n:
                    ax.axvline(offset - 2, color="grey", linestyle="--", linewidth=0.5)
                variable = rapidity_labels[bins["rapidity_variable"].iloc[0]]
                ax.text(
                    x.mean(),
                    0.96,
                    f"${variable}={rapidity:g}$",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=9,
                )
                take = np.arange(0, len(bins), max(1, int(np.ceil(len(bins) / 4))))
                ticks.extend(x[take])
                labels.extend(f"{value:g}" for value in bins["scale"].iloc[take])
                offset += len(bins) + 3
            ax.axhline(0, color="grey", linewidth=0.5)
            ax.set_xticks(ticks, labels, rotation=45, ha="right", fontsize=8)
            scale = r"$p_T$" if row["scale_variable"].iloc[0] == "pT" else r"$m_{jj}$"
            ax.set(xlabel=scale + " [GeV]", ylabel="Correction / relative error [%]")
            if pd.notna(yb):
                ax.set_title(f"$y_b={yb:g}$")
            ax.margins(x=0.01, y=0.2)
            ax.legend(loc="upper right", bbox_to_anchor=(1, -0.22), ncol=3, fontsize=8)
        fig.tight_layout()
        yield fig


pc_fit_jet_correction_comparison = collect("pc_jet_correction_comparison", ("pc_fit_context",))


@table
def pc_fits_jet_correction_comparison(fits, pc_merge_api):
    """Resolve fits sequentially, retaining only compact pointwise comparison tables."""
    check(len(fits) > 0, "Provide at least one fit for the jet correction comparison.")
    check(len({fit.name for fit in fits}) == len(fits), "Comparison fits must be distinct.")
    return pd.concat(
        {fit.name: pc_merge_api.pc_fit_jet_correction_comparison(fit=fit.name)[0] for fit in fits},
        names=["fit"],
    )


@table
def pc_fits_jet_correction_exceedance(pc_fits_jet_correction_comparison):
    """Per-fit rapidity-bin summaries of the correction-size criterion."""
    if pc_fits_jet_correction_comparison.empty:
        empty = pc_jet_correction_exceedance(pc_fits_jet_correction_comparison.droplevel("fit"))
        return pd.concat({"": empty}, names=["fit"])
    return pd.concat(
        {
            name: pc_jet_correction_exceedance(points.droplevel("fit"))
            for name, points in pc_fits_jet_correction_comparison.groupby(level="fit", sort=False)
        },
        names=["fit"],
    )


@figuregen
def plot_pc_fits_jet_correction_comparison(pc_fits_jet_correction_comparison):
    """One dataset figure per fit, using each fit's own posterior and cuts."""
    for name, points in pc_fits_jet_correction_comparison.groupby(level="fit", sort=False):
        for fig in plot_pc_jet_correction_comparison(points.droplevel("fit")):
            fig.suptitle(f"{name}\n{fig.get_suptitle()}")
            fig.tight_layout()
            yield fig
