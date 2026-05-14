#!/usr/bin/env python3
"""Run NNPDF Shapley value analysis from a YAML runcard.

Usage
-----
    vp-shapley runcards/sv_dis_hera.yaml
    vp-shapley runcards/sv_global.yaml --output results/sv
"""

import sys
sys.path.insert(0, '/home/daksh/nnpdf/nnpdf/validphys2/src')

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
import yaml

import lhapdf
from validphys.shapley.setup import setup_observables
from validphys.shapley.analyzer import NNPDFShapleyAnalyzer
from validphys.shapley.perturbation import apply_gaussian_perturbation
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

from shapley_values import save_results
from validphys.shapley import perturbation

_MUTED_SERIES_COLORS = [
    "#E45756",
    "#F58518",
    "#EECA3B",
    "#71BD00B5",
    "#72B7B2",
    "#4C78A8",
    "#B279A2",
    "#6F4C9B",
    "#BAB0AC"
]

_SERIES_MARKERS = ["o", "s", "D", "^", "v", "P", "X", "<", ">", "h"]
_SOFT_POSITIVE_COLOR = "#518500"
_SOFT_NEGATIVE_COLOR = "#FFBF00"


def _plot_sv_bar(sv, labels, title=None, sv_err=None,
                 ax=None, figsize=(12, 6),
                 ylabel="Shapley Value (delta chi2/N)",
                 positive_color=_SOFT_POSITIVE_COLOR,
                 negative_color=_SOFT_NEGATIVE_COLOR,
                 alpha=0.85, show_legend=True):
    """Bar chart of Shapley values with optional symmetric error bars.

    Parameters
    ----------
    sv : array-like, shape (n,)
    labels : list of str
    title : str, optional
    sv_err : array-like, shape (n,) or None
        1-sigma uncertainty per flavour (e.g. replica std or std-error).
        When provided, symmetric black error bars are drawn on each bar.
    ax : matplotlib Axes, optional
        Draw into this axes; otherwise a new figure is created.
    figsize, ylabel, positive_color, negative_color, alpha
        Forwarded to the underlying bar plot.
    """
    sv = np.asarray(sv, dtype=float)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    bars = ax.bar(labels, sv)
    for bar, val in zip(bars, sv):
        bar.set_color(positive_color if val >= 0 else negative_color)
        bar.set_alpha(alpha)

    if sv_err is not None:
        sv_err = np.asarray(sv_err, dtype=float)
        if np.any(np.isfinite(sv_err)):
            sv_err = np.where(np.isfinite(sv_err), sv_err, 0.0)
            ax.errorbar(
                labels, sv,
                yerr=sv_err,
                fmt="none",
                ecolor="#333333",
                elinewidth=1.1,
                capsize=3.2,
                capthick=1.1,
                zorder=5,
            )

    ax.axhline(0, color="#444444", lw=0.8, alpha=0.9)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", ls="--", alpha=0.35, linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    if title:
        ax.set_title(title)

    handles = [
        mpatches.Patch(color=positive_color, alpha=alpha, label="SV > 0"),
        mpatches.Patch(color=negative_color, alpha=alpha, label="SV < 0"),
    ]
    if sv_err is not None and np.any(np.isfinite(sv_err)):
        handles.append(
            mlines.Line2D([], [], color="#333333", linewidth=1.1, label="1\u03c3 std")
        )
    if show_legend:
        ax.legend(handles=handles, loc="upper right", frameon=False)
    fig.tight_layout()
    return fig


def load_runcard(path):
    """Load and validate a YAML runcard."""
    with open(path) as f:
        cfg = yaml.safe_load(f)

    required = ["pdf_name", "datasets", "experiments"]
    for key in required:
        if key not in cfg:
            raise KeyError(f"Missing required key '{key}' in runcard")

    if not isinstance(cfg["experiments"], (list, dict)):
        raise TypeError("'experiments' must be a list or a mapping")

    return cfg

def compute_ceilings(cfg, x_min = 1e-5, x_max = 0.9999, 
                    x_points = 500, scaling='log', nf = 6, ubound=50):
    """
    Compute and add a ceilings list to each experiment in the cfg file

    Parameters:
    -----------
    cfg: dict
        Configuration file from the runcard
    x_min: float
        Min. x of PDF evaluation
    x_max: float 
        Max. x of PDF evaluation
    x_points: int
        Total points of PDF evaluation
    scaling: str
        Scaling for x values for PDF evaluation
    nf: int
        No. of quark flavors 
    
    Returns:
    --------
        cfg: dict
            Returns cfg file with amplitude ceilings added to each experiment
    """

    # Build PDF
    pdf_dict = _build_setup_context(cfg)
    pdf = lhapdf.mkPDF(str(pdf_dict['pdf']), 0)
    Q0 = pdf_dict['observables'][0].Q0

    # Compute amp ceiling
    pids = [1,2,3,4,5,6]

    gluon_ceiling_plus, quark_ceiling_plus = float('inf'), float('inf')
    gluon_ceiling_minus, quark_ceiling_minus = float('inf'), float('inf')*np.ones(nf) 
    
    xgrid = np.logspace(np.log10(x_min), np.log10(x_max), x_points)
    if scaling == 'log':
        def pert_exponent(x,mu,sig): return 0.5 * ((np.log10(x)-np.log10(mu))/sig) ** 2
    elif scaling == 'linear':
        def pert_exponent(x,mu,sig): return 0.5 * ((x-mu)/sig) ** 2
    else:
        raise ValueError(f"Unknown scaling '{scaling}'. Scaling must be either linear or logarithmic.")
    
    all_exp = cfg['experiments']
    for exp in all_exp:
        pert = exp['perturbation']
        mu = pert['mu']
        sig = pert['sigma']
        mode = pert['mode']
        if mode == 'physical':
            product = np.zeros(nf+2) # nf quark flavor, 1 gluon, 1 singlet
            for x in xgrid:
                # Get x*f(x, Q^2) 
                xf_singlet = sum(pdf.xfxQ(i, x, Q0) + pdf.xfxQ(-i, x, Q0) for i in pids[:nf])
                xg = pdf.xfxQ(21,x,Q0) # gluon pid = 21
                xquark = [pdf.xfxQ(i,x,Q0) for i in pids[:nf]]
                
                # Compute product and track minimums
                exponent = pert_exponent(x,mu,sig)
                exponent = np.clip(exponent, 0, ubound)   # can clip at 50 as min. x*f(x) goes to order -7 to -10 for Q=1.65 to 200GeV
                gaussian = np.exp(exponent)

                # Only consider physical PDF values (>0)
                # Otherwise the ceiling turns out to be negative
                xPDFs = [xf_singlet] + [xg] + xquark
                
                for idx, xPDF in enumerate(xPDFs):
                    if xPDF > 0:
                        product[idx] = xPDF * gaussian

                xsingletG = product[0]
                xgluonG = 0.5 * product[1]
                xquarkG = product[1:]

                gluon_ceiling_plus = min(gluon_ceiling_plus, xsingletG)
                quark_ceiling_plus = min(quark_ceiling_plus, xgluonG)
                quark_ceiling_minus = [min(quark_ceiling_minus[i-1], xquarkG[i-1]) for i in pids[:nf]]
            
            gluon_ceiling_minus = 2 * quark_ceiling_plus

            # To make life easier, we save all ceilings in the same format as the dict 'flavor_info'
            ceilings_plus = quark_ceiling_plus*np.ones(2*nf+1)
            ceilings_plus[nf] = gluon_ceiling_plus   # since gluon index in flavor_info is 6-(6-nf)
            ceilings_minus = quark_ceiling_minus[::-1] + [gluon_ceiling_minus] + quark_ceiling_minus
        else:
            return cfg
        # Now add ceilings to each perturbation
        pert['ceilings_plus'] = ceilings_plus
        pert['ceilings_minus'] = ceilings_minus

    return cfg

def _safe_token(value):
    """Build filesystem-safe tokens for output naming."""
    token = str(value).strip().replace(".", "p")
    token = re.sub(r"[^A-Za-z0-9_-]+", "_", token)
    return token.strip("_")


def _resolve_output_dir(cfg, runcard_path, cli_output):
    """Resolve output directory with precedence: CLI > runcard > default."""
    shapley_root = Path(__file__).resolve().parent.parent
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    if cli_output:
        out = Path(cli_output)
        return out if out.is_absolute() else (shapley_root / out)

    if cfg.get("output_dir"):
        out = Path(cfg["output_dir"])
        return out if out.is_absolute() else (shapley_root / out)

    output_root_cfg = cfg.get("output_root", str(shapley_root / "sv_results"))
    output_root = Path(output_root_cfg)
    if not output_root.is_absolute():
        output_root = shapley_root / output_root

    run_name = cfg.get("run_name")
    label = _safe_token(run_name) if run_name else _safe_token(Path(runcard_path).stem)
    return output_root / f"{ts}_{label}"


def _normalize_experiments(cfg):
    """Normalize experiments into a list of (name, config) tuples."""
    experiments = cfg["experiments"]
    if isinstance(experiments, dict):
        items = []
        for name, exp_cfg in experiments.items():
            merged = dict(exp_cfg or {})
            merged["name"] = merged.get("name", name)
            items.append(merged)
    else:
        items = experiments

    normalized = []
    for i, exp in enumerate(items, start=1):
        if not isinstance(exp, dict):
            raise TypeError(f"Experiment #{i} must be a mapping")

        exp_name = str(exp.get("name", f"exp_{i:02d}"))
        pert = exp.get("perturbation")
        if pert is None:
            raise KeyError(f"Missing 'perturbation' in experiment '{exp_name}'")
        for key in ["mu", "sigma", "amplitude"]:
            if key not in pert:
                raise KeyError(f"Missing {key} in perturbation for experiment '{exp_name}'")

        exp_cfg = dict(cfg)
        exp_cfg["perturbation"] = pert
        for key in [
            "basis",
            "n_jobs",
            "enforce_sumrules",
            "n_replicas",
            "n_sign_samples",
            "random_seed",
            "run_name",
            "stabilization",
        ]:
            if key in exp:
                exp_cfg[key] = exp[key]

        normalized.append((exp_name, exp_cfg))

    return normalized


def _build_setup_context(cfg):
    """Load PDF and observables once for all experiments in a runcard."""
    pdf, observables, flavor_info, _ = setup_observables(
        pdf_name=cfg["pdf_name"],
        datasets=cfg["datasets"],
        theory_id=cfg.get("theory_id", 708),
        use_cuts=cfg.get("use_cuts", "internal"),
        variant=cfg.get("variant", None),
    )
    return {
        "pdf": pdf,
        "observables": observables,
        "flavor_info": flavor_info,
    }


def _format_axis_tick(val):
    """Compact numeric tick labels suitable for x-scan values."""
    val = float(val)
    if val == 0.0:
        return "0"
    abs_v = abs(val)
    if 1e-2 <= abs_v < 1e2:
        return f"{val:g}"
    text = f"{val:.0e}"
    return text.replace("e-0", "e-").replace("e+0", "e+")


def _resolve_scan_axis(exp_names, experiment_meta, evenly_spaced_bins=True):
    """Return axis positions and display labels for experiment comparisons."""
    mu_vals = []
    xspaces = []
    for exp_name in exp_names:
        pert = (experiment_meta.get(exp_name) or {}).get("perturbation", {})
        try:
            mu_vals.append(float(pert.get("mu")))
        except (TypeError, ValueError):
            mu_vals = None
            break
        xspaces.append(str(pert.get("xspace", "linear")).strip().lower())

    n_exp = len(exp_names)
    if not mu_vals or len(mu_vals) != n_exp:
        x = np.arange(n_exp, dtype=float)
        return {
            "order": np.arange(n_exp),
            "x": x,
            "tick_labels": exp_names,
            "xscale": "linear",
            "xlabel": "Experiment",
            "is_index_axis": True,
        }

    order = np.argsort(mu_vals)
    x_sorted = np.asarray(mu_vals, dtype=float)[order]
    if evenly_spaced_bins:
        # Use evenly spaced bins for readability while keeping sorted mu labels.
        x_bins = np.arange(n_exp, dtype=float)
        return {
            "order": order,
            "x": x_bins,
            "tick_labels": [_format_axis_tick(v) for v in x_sorted],
            "xscale": "linear",
            "xlabel": "Perturbation center x (mu, evenly spaced bins)",
            "is_index_axis": True,
        }

    all_positive = bool(np.all(x_sorted > 0))
    use_log = all_positive and all(xs == "logx" for xs in xspaces)
    return {
        "order": order,
        "x": x_sorted,
        "tick_labels": [_format_axis_tick(v) for v in x_sorted],
        "xscale": "log" if use_log else "linear",
        "xlabel": "Perturbation center x (mu)",
        "is_index_axis": False,
    }


def _plot_sv_scan_comparison(
    sv_matrix,
    err_matrix,
    labels,
    basis,
    output_dir,
    axis_info,
    output_stem,
    title_suffix=None,
    dodge_step_pt=4.5,
    use_dodge=True,
):
    """Plot SV vs scan-x with per-flavor markers and uncertainty bars."""
    x = np.asarray(axis_info["x"], dtype=float)
    use_index_axis = bool(axis_info["is_index_axis"])

    with matplotlib.rc_context(
        {
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 14,
            "legend.fontsize": 10,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    ):
        fig, ax = plt.subplots(figsize=(10.5, 6.2))
        n_series = len(labels)
        if (not use_dodge) or n_series <= 1:
            dodge_offsets_pt = np.zeros(n_series, dtype=float)
        else:
            # Constant display-space dodge keeps interpretation of x intact
            # for both linear and logarithmic axes.
            half_span_pt = 0.5 * dodge_step_pt * (n_series - 1)
            dodge_offsets_pt = np.linspace(
                -half_span_pt, half_span_pt, n_series, dtype=float
            )

        for i, flav in enumerate(labels):
            y = sv_matrix[:, i]
            yerr = err_matrix[:, i] if err_matrix is not None else None
            valid = np.isfinite(y)
            if not np.any(valid):
                continue

            xv = x[valid]
            yv = y[valid]
            color = _MUTED_SERIES_COLORS[i % len(_MUTED_SERIES_COLORS)]
            marker = _SERIES_MARKERS[i % len(_SERIES_MARKERS)]

            if yerr is not None:
                yerrv = np.asarray(yerr[valid], dtype=float)
                yerrv = np.where(np.isfinite(yerrv), yerrv, 0.0)
            else:
                yerrv = None

            err_container = ax.errorbar(
                xv,
                yv,
                yerr=yerrv,
                linestyle="none",
                linewidth=0.0,
                marker=marker,
                markersize=5.0,
                markerfacecolor="white",
                markeredgewidth=1.0,
                capsize=2.8,
                elinewidth=1.0,
                color=color,
                alpha=0.95,
                label=flav,
                zorder=3,
            )

            dx_pt = float(dodge_offsets_pt[i])
            if dx_pt != 0.0:
                shift = mtransforms.ScaledTranslation(
                    dx_pt / 72.0, 0.0, fig.dpi_scale_trans
                )
                shifted_data = ax.transData + shift
                data_line, cap_lines, bar_linecols = err_container.lines
                if data_line is not None:
                    data_line.set_transform(shifted_data)
                for artist in cap_lines:
                    artist.set_transform(shifted_data)
                for artist in bar_linecols:
                    artist.set_transform(shifted_data)

        if axis_info["xscale"] == "log":
            ax.set_xscale("log")

        ax.set_xticks(x)
        ax.set_xticklabels(axis_info["tick_labels"], rotation=0 if not use_index_axis else 35)
        ax.set_xlabel(axis_info["xlabel"])
        ax.set_ylabel("Shapley Value (delta chi2/N)")
        title = f"Shapley comparison ({basis} basis)"
        if title_suffix:
            title = f"{title} - {title_suffix}"
        ax.set_title(title)
        ax.axhline(0.0, color="#555555", linewidth=0.9, alpha=0.9)
        ax.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.35)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=False,
            ncol=1,
            borderaxespad=0.0,
        )

        fig.tight_layout(rect=[0, 0, 0.82, 1])

        png_path = output_dir / f"{output_stem}_{basis}.png"
        pdf_path = output_dir / f"{output_stem}_{basis}.pdf"
        fig.savefig(png_path, dpi=220, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved comparison plot: {png_path}")
        print(f"Saved comparison plot: {pdf_path}")


def _format_experiment_panel_title(exp_name, exp_meta):
    """Compact per-panel title for bar comparison figures."""
    pert = (exp_meta or {}).get("perturbation", {})
    mu = pert.get("mu")
    sigma = pert.get("sigma")
    if mu is None:
        return str(exp_name)
    mu_txt = _format_axis_tick(mu)
    if sigma is None:
        return f"{exp_name} | mu={mu_txt}"
    return f"{exp_name} | mu={mu_txt}, sigma={sigma}"


def _plot_sv_bar_comparison(
    sv_matrix,
    err_matrix,
    labels,
    basis,
    exp_names_ordered,
    output_dir,
    experiment_meta,
    output_stem="shapley_comparison_bars",
    title_tag="even",
):
    """Multi-panel bar comparisons across experiments."""
    n = len(exp_names_ordered)
    if n == 0:
        return

    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))

    with matplotlib.rc_context(
        {
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    ):
        fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 4.4 * nrows))
        axes_flat = np.atleast_1d(axes).ravel()

        for i, (ax, exp_name) in enumerate(zip(axes_flat, exp_names_ordered)):
            title = _format_experiment_panel_title(
                exp_name, experiment_meta.get(exp_name, {})
            )
            _plot_sv_bar(
                sv_matrix[i, :],
                labels,
                title=title,
                sv_err=err_matrix[i, :],
                ax=ax,
                show_legend=False,
            )
            if (i % ncols) != 0:
                ax.set_ylabel("")

        for ax in axes_flat[n:]:
            ax.set_visible(False)

        fig.suptitle(
            f"Shapley bar comparison ({basis} basis, {title_tag})", y=0.995
        )

        fig.tight_layout(rect=[0, 0, 1, 0.98])

        png_path = output_dir / f"{output_stem}_{basis}.png"
        pdf_path = output_dir / f"{output_stem}_{basis}.pdf"
        fig.savefig(png_path, dpi=220, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved comparison plot: {png_path}")
        print(f"Saved comparison plot: {pdf_path}")


def _save_comparison_plots(all_experiment_results, output_dir, experiment_meta):
    """Save per-basis SV-vs-x and bar-comparison plots across experiments."""
    if len(all_experiment_results) < 2:
        return

    basis_names = sorted(
        {
            basis
            for exp_results in all_experiment_results.values()
            for basis in exp_results.keys()
        }
    )
    exp_names = list(all_experiment_results.keys())
    axis_info_bins = _resolve_scan_axis(
        exp_names, experiment_meta, evenly_spaced_bins=True
    )
    axis_info_truex = _resolve_scan_axis(
        exp_names, experiment_meta, evenly_spaced_bins=False
    )
    order_bins = axis_info_bins["order"]
    exp_names_ordered = [exp_names[i] for i in order_bins]

    for basis in basis_names:
        basis_entries = [
            all_experiment_results[exp_name].get(basis, {})
            for exp_name in exp_names
            if basis in all_experiment_results[exp_name]
        ]
        deterministic_basis = bool(basis_entries) and all(
            bool(entry.get("deterministic_sign_symmetrized", False))
            for entry in basis_entries
        )
        replica_basis = bool(basis_entries) and all(
            bool(entry.get("per_replica", False))
            for entry in basis_entries
        )
        if deterministic_basis and replica_basis:
            mode_tag = "det_rep"
        elif deterministic_basis:
            mode_tag = "det"
        elif replica_basis:
            mode_tag = "rep"
        else:
            mode_tag = "pert"

        def _plot_game(value_key, std_key, err_key, stem_tag, title_tag):
            labels = None
            for exp_name in exp_names:
                if basis in all_experiment_results[exp_name]:
                    labels = list(
                        (all_experiment_results[exp_name][basis].get(value_key) or {}).keys()
                    )
                    if labels:
                        break
            if labels is None:
                return

            sv_rows = []
            err_rows = []
            for exp_name in exp_names:
                basis_results = all_experiment_results[exp_name].get(basis, {})
                sv_map = basis_results.get(value_key, {}) or {}
                err_map = (
                    basis_results.get(std_key)
                    or basis_results.get(err_key)
                    or {}
                )
                sv_rows.append(
                    np.array([float(sv_map.get(lbl, 0.0)) for lbl in labels])
                )
                err_rows.append(
                    np.array([float(err_map.get(lbl, np.nan)) for lbl in labels])
                    if err_map else np.full(len(labels), np.nan)
                )

            sv_matrix_bins = np.asarray(sv_rows, dtype=float)[order_bins, :]
            err_matrix_bins = np.asarray(err_rows, dtype=float)[order_bins, :]
            _plot_sv_scan_comparison(
                sv_matrix=sv_matrix_bins,
                err_matrix=err_matrix_bins,
                labels=labels,
                basis=basis,
                output_dir=output_dir,
                axis_info=axis_info_bins,
                output_stem=f"shapley_comparison_{stem_tag}",
                title_suffix=f"{title_tag}, binned x, dodged points",
                dodge_step_pt=4.5,
                use_dodge=True,
            )

            order_truex = axis_info_truex["order"]
            sv_matrix_truex = np.asarray(sv_rows, dtype=float)[order_truex, :]
            err_matrix_truex = np.asarray(err_rows, dtype=float)[order_truex, :]
            _plot_sv_scan_comparison(
                sv_matrix=sv_matrix_truex,
                err_matrix=err_matrix_truex,
                labels=labels,
                basis=basis,
                output_dir=output_dir,
                axis_info=axis_info_truex,
                output_stem=f"shapley_comparison_truex_{stem_tag}",
                title_suffix=f"{title_tag}, true x, overlapping",
                use_dodge=False,
            )

            _plot_sv_bar_comparison(
                sv_matrix=sv_matrix_bins,
                err_matrix=err_matrix_bins,
                labels=labels,
                basis=basis,
                exp_names_ordered=exp_names_ordered,
                output_dir=output_dir,
                experiment_meta=experiment_meta,
                output_stem=f"shapley_comparison_bars_{stem_tag}",
                title_tag=title_tag,
            )

        _plot_game(
            value_key="shapley_values",
            std_key="shapley_std",
            err_key="shapley_err",
            stem_tag=mode_tag,
            title_tag=mode_tag,
        )

        _save_dense_flavor_heatmap(
            all_experiment_results=all_experiment_results,
            output_dir=output_dir,
            experiment_meta=experiment_meta,
            basis=basis,
            mode_tag=mode_tag,
        )


def _compute_log_bin_edges(xvals):
    """Compute log-spaced bin edges around sorted positive x centers."""
    x = np.asarray(xvals, dtype=float)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("xvals must be a non-empty 1D array")
    if not np.all(x > 0.0):
        raise ValueError("xvals must be strictly positive for log edges")

    if x.size == 1:
        return np.asarray([x[0] / np.sqrt(10.0), x[0] * np.sqrt(10.0)], dtype=float)

    logx = np.log10(x)
    mids = 0.5 * (logx[:-1] + logx[1:])
    edges = np.empty(x.size + 1, dtype=float)
    edges[1:-1] = mids
    edges[0] = logx[0] - (mids[0] - logx[0])
    edges[-1] = logx[-1] + (logx[-1] - mids[-1])
    return np.power(10.0, edges)


def _save_dense_flavor_heatmap(
    all_experiment_results,
    output_dir,
    experiment_meta,
    basis,
    mode_tag="pert",
):
    """Save dense flavor x-scan heatmap (+ global |phi| band)."""
    # Keep this output focused on dense scans so existing coarse workflows stay unchanged.
    if len(all_experiment_results) < 10:
        return

    exp_names = list(all_experiment_results.keys())
    usable = []
    for exp_name in exp_names:
        basis_results = all_experiment_results.get(exp_name, {}).get(basis)
        if not basis_results:
            continue
        pert = (experiment_meta.get(exp_name) or {}).get("perturbation", {})
        try:
            mu = float(pert.get("mu"))
        except (TypeError, ValueError):
            continue

        phi_even_map = basis_results.get("phi_even") or basis_results.get("shapley_values") or {}
        if not phi_even_map:
            continue
        usable.append((exp_name, mu, phi_even_map))

    if len(usable) < 2:
        return

    usable.sort(key=lambda item: item[1])
    mu_sorted = np.asarray([item[1] for item in usable], dtype=float)
    labels = list(usable[0][2].keys())
    n_labels = len(labels)

    phi_matrix = np.full((n_labels, len(usable)), np.nan, dtype=float)
    for col, (_, _, phi_map) in enumerate(usable):
        for row, label in enumerate(labels):
            val = phi_map.get(label)
            if val is not None:
                phi_matrix[row, col] = float(val)

    global_importance = np.nanmean(np.abs(phi_matrix), axis=0)
    if not np.any(np.isfinite(phi_matrix)):
        return

    x_edges = _compute_log_bin_edges(mu_sorted)
    y_edges_top = np.arange(n_labels + 1, dtype=float)
    y_edges_bottom = np.array([0.0, 1.0], dtype=float)

    with matplotlib.rc_context(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    ):
        fig, (ax_top, ax_bottom) = plt.subplots(
            2,
            1,
            figsize=(12.5, 6.6),
            sharex=True,
            gridspec_kw={"height_ratios": [9.0, 1.35], "hspace": 0.08},
        )

        finite_abs = np.abs(phi_matrix[np.isfinite(phi_matrix)])
        absmax = float(np.max(finite_abs)) if finite_abs.size else 0.0
        if absmax <= 0.0:
            absmax = 1e-12
        nonzero_abs = finite_abs[finite_abs > 0.0]
        if nonzero_abs.size:
            # Symlog keeps sign while expanding contrast around small magnitudes.
            linthresh = float(np.percentile(nonzero_abs, 25.0))
            linthresh = max(min(linthresh, absmax), 1e-12)
        else:
            linthresh = 1e-12

        top_mesh = ax_top.pcolormesh(
            x_edges,
            y_edges_top,
            phi_matrix,
            cmap="RdBu_r",
            norm=SymLogNorm(
                linthresh=linthresh,
                linscale=1.0,
                vmin=-absmax,
                vmax=absmax,
                base=10.0,
            ),
            shading="auto",
        )

        global_2d = global_importance[np.newaxis, :]
        gmax = float(np.nanmax(global_2d)) if np.any(np.isfinite(global_2d)) else 0.0
        if gmax <= 0.0:
            gmax = 1e-12
        pos_global = global_2d[np.isfinite(global_2d) & (global_2d > 0.0)]
        gmin = float(np.min(pos_global)) if pos_global.size else 1e-12
        gmin = min(gmin, gmax)
        bot_mesh = ax_bottom.pcolormesh(
            x_edges,
            y_edges_bottom,
            global_2d,
            cmap="YlGnBu",
            norm=LogNorm(vmin=max(gmin, 1e-12), vmax=max(gmax, gmin, 1e-12)),
            shading="auto",
        )

        ax_top.set_xscale("log")
        ax_bottom.set_xscale("log")
        ax_top.set_xlim(float(x_edges[0]), float(x_edges[-1]))

        ax_top.set_yticks(np.arange(n_labels, dtype=float) + 0.5)
        ax_top.set_yticklabels(labels)
        ax_top.invert_yaxis()
        ax_top.set_ylabel("flavor")
        top_tag = f"signed {mode_tag}"
        ax_top.set_title(f"Dense x-scan heatmap ({top_tag})")

        ax_bottom.set_yticks([0.5])
        ax_bottom.set_yticklabels(["global"])
        ax_bottom.invert_yaxis()
        ax_bottom.set_ylabel("band")

        # Use true mu geometry (log-scaled x edges). To keep labels readable,
        # only show a subset of exact mu tick labels.
        n_ticks = min(10, len(mu_sorted))
        tick_idx = np.unique(np.linspace(0, len(mu_sorted) - 1, num=n_ticks, dtype=int))
        tick_vals = mu_sorted[tick_idx]
        ax_bottom.set_xticks(tick_vals)
        ax_bottom.set_xticklabels([_format_axis_tick(v) for v in tick_vals], rotation=30, ha="right")
        ax_bottom.set_xlabel("perturbation center x (mu)")

        cbar_top = fig.colorbar(top_mesh, ax=ax_top, pad=0.015)
        cbar_top.set_label(f"signed {mode_tag} (symlog)")
        cbar_bottom = fig.colorbar(bot_mesh, ax=ax_bottom, pad=0.015)
        cbar_bottom.set_label("mean_j |phi_j(x)| (log)")

        fig.subplots_adjust(left=0.08, right=0.92, bottom=0.14, top=0.94, hspace=0.08)

        stem = f"shapley_dense_heatmap_{mode_tag}"
        png_path = output_dir / f"{stem}_{basis}.png"
        pdf_path = output_dir / f"{stem}_{basis}.pdf"
        fig.savefig(png_path, dpi=220, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved comparison plot: {png_path}")
        print(f"Saved comparison plot: {pdf_path}")


def _save_diagnostic_files(diag, output_dir, basis):
    """Write per-coalition chi2 CSV, marginal-contribution CSV, and stats JSON.

    Files produced
    --------------
    coalition_chi2_<basis>.csv
        One row per coalition sorted by chi2 (descending).
        Columns: coalition_labels, size, chi2, is_outlier

    marginal_contributions_<basis>.csv
        One row per (player, coalition-without-player) pair.
        Columns: player, coalition_without, v_without, v_with, delta_v, is_outlier

    diagnostic_stats_<basis>.json
        Summary statistics for chi2 and per-player marginals.
    """
    import copy

    # -- coalition_chi2 CSV -------------------------------------------------
    outlier_threshold = diag.get("outlier_chi2_threshold", float("inf"))
    outlier_set = {
        tuple(oc["coalition"]) for oc in diag.get("outlier_coalitions", [])
    }
    # Rebuild from the marginal list to get all coalitions (including empty).
    # Also collect the full coalition set from chi2_stats (n_coalitions count).
    # We rely on _marginal_contributions for the data.
    marginals = diag.get("_marginal_contributions", [])

    # Build dict coalition -> chi2 from marginals (v_without is chi2 of coalition).
    chi2_from_marginals = {}
    for m in marginals:
        coal = tuple(m["coalition_without"])
        chi2_from_marginals[coal] = m["v_without"]
        coal_with = tuple(sorted(m["coalition_without"] + [m["player_idx"]]))
        chi2_from_marginals[coal_with] = m["v_with"]

    # Sort by chi2 descending.
    sorted_coalitions = sorted(
        chi2_from_marginals.items(), key=lambda kv: -kv[1]
    )

    # Map player indices -> short labels (for human-readable coalition column).
    # We need the flavor_short list; it isn't in diag, so re-derive from marginals.
    idx_to_label = {}
    for m in marginals:
        idx_to_label[m["player_idx"]] = m["player"]

    coalition_csv_path = output_dir / f"coalition_chi2_{basis}.csv"
    with open(coalition_csv_path, "w") as f:
        f.write("coalition_labels,size,chi2,is_outlier\n")
        for coal, chi2 in sorted_coalitions:
            labels_str = "|".join(idx_to_label.get(i, str(i)) for i in coal)
            is_out = 1 if coal in outlier_set or chi2 > outlier_threshold else 0
            f.write(f"{labels_str},{len(coal)},{chi2:.8f},{is_out}\n")
    print(f"Saved: {coalition_csv_path}")

    # -- marginal_contributions CSV -----------------------------------------
    contrib_csv_path = output_dir / f"marginal_contributions_{basis}.csv"
    with open(contrib_csv_path, "w") as f:
        f.write("player,coalition_without,v_without,v_with,delta_v,is_outlier\n")
        for m in sorted(marginals, key=lambda x: -abs(x["delta_v"])):
            coal_str = "|".join(
                idx_to_label.get(i, str(i)) for i in m["coalition_without"]
            )
            f.write(
                f"{m['player']},{coal_str},"
                f"{m['v_without']:.8f},{m['v_with']:.8f},"
                f"{m['delta_v']:.8f},{m.get('is_outlier', 0)}\n"
            )
    print(f"Saved: {contrib_csv_path}")

    # -- diagnostic_stats JSON ----------------------------------------------
    # Strip the internal key before serialising.
    diag_public = {k: v for k, v in diag.items() if k != "_marginal_contributions"}
    stats_path = output_dir / f"diagnostic_stats_{basis}.json"
    with open(stats_path, "w") as f:
        json.dump(diag_public, f, indent=2)
    print(f"Saved: {stats_path}")


def _write_shapley_csv(path, labels, values):
    """Write Shapley values to CSV with standard schema."""
    with open(path, "w") as f:
        f.write("flavour,shapley_value\n")
        for lbl, val in zip(labels, values):
            f.write(f"{lbl},{float(val):.8f}\n")


def _write_shapley_per_replica_csv(path, labels, values_per_replica):
    """Write per-replica Shapley values as a replica-by-flavour matrix."""
    matrix = np.asarray(values_per_replica, dtype=float)
    with open(path, "w") as f:
        f.write("replica," + ",".join(labels) + "\n")
        for replica_idx, row in enumerate(matrix, start=1):
            values_txt = ",".join(f"{float(val):.8f}" for val in row)
            f.write(f"{replica_idx},{values_txt}\n")


def _write_shapley_per_sample_csv(path, labels, values_per_sample):
    """Write per-sign-sample Shapley values as a sample-by-flavour matrix."""
    matrix = np.asarray(values_per_sample, dtype=float)
    with open(path, "w") as f:
        f.write("sample," + ",".join(labels) + "\n")
        for sample_idx, row in enumerate(matrix, start=1):
            values_txt = ",".join(f"{float(val):.8f}" for val in row)
            f.write(f"{sample_idx},{values_txt}\n")


def _serialize_shapley_per_replica(labels, values_per_replica):
    """Convert per-replica SV matrix into JSON-native flavour -> list form."""
    if values_per_replica is None:
        return None
    matrix = np.asarray(values_per_replica, dtype=float)
    return {
        lbl: [float(val) for val in matrix[:, i]]
        for i, lbl in enumerate(labels)
    }


def _serialize_shapley_per_sample(labels, values_per_sample):
    """Convert per-sign-sample SV matrix into JSON-native flavour -> list form."""
    if values_per_sample is None:
        return None
    matrix = np.asarray(values_per_sample, dtype=float)
    return {
        lbl: [float(val) for val in matrix[:, i]]
        for i, lbl in enumerate(labels)
    }


def _resolve_stabilization_cfg(cfg):
    """Resolve stabilization options with safe defaults."""
    stab = cfg.get("stabilization", {}) or {}
    if not isinstance(stab, dict):
        raise TypeError("'stabilization' must be a mapping if provided")

    action = str(stab.get("action", "exclude_dataset")).strip().lower()
    if action not in {"exclude_dataset", "report_only"}:
        raise ValueError(
            "stabilization.action must be 'exclude_dataset' or 'report_only'"
        )

    threshold = float(stab.get("dataset_delta_chi2_threshold", 1e4))
    if threshold <= 0:
        raise ValueError("stabilization.dataset_delta_chi2_threshold must be > 0")

    max_outlier_coalitions = int(stab.get("max_outlier_coalitions", 5))
    if max_outlier_coalitions < 1:
        raise ValueError("stabilization.max_outlier_coalitions must be >= 1")

    return {
        # Native behavior unless explicitly disabled in the runcard.
        "enabled": bool(stab.get("enabled", True)),
        "action": action,
        "dataset_delta_chi2_threshold": threshold,
        "max_outlier_coalitions": max_outlier_coalitions,
        "rerun_stable": bool(stab.get("rerun_stable", True)),
    }


def _dataset_mean_chi2_for_coalition(analyzer, coalition, pert, sign_matrices=None):
    """Compute per-dataset mean chi2 for a single coalition."""
    mu = pert["mu"]
    sigma = pert["sigma"]
    amplitude = pert["amplitude"]
    mode = pert.get("mode", "additive")
    xspace = pert.get("xspace", "linear")
    sign_matrices = [None] if sign_matrices is None else list(sign_matrices)

    rows = []
    for obs in analyzer.observables:
        mean_samples = []
        for sign_matrix in sign_matrices:
            sr_norm = None
            if analyzer.enforce_sumrules:
                sr_norm = analyzer._compute_sumrule_norm(
                    coalition, mu, sigma, amplitude, mode, xspace,
                    random_sign_matrix=sign_matrix,
                )

            if analyzer.basis == "flavor":
                gv_pert_list = []
                perturb_idx = [analyzer.flavor_indices[p] for p in coalition]
                perturb_signs = (
                    sign_matrix[:, coalition] if sign_matrix is not None else None
                )
                for idx, entry in enumerate(obs.fk_entries):
                    gv_flav = analyzer._get_flavor_gv_for_entry(obs, idx)
                    gv_flav_calib = analyzer._get_calibration_flavor_gv_for_entry(obs, idx)
                    gv_pert = apply_gaussian_perturbation(
                        gv_flav,
                        perturb_idx,
                        mu,
                        sigma,
                        amplitude,
                        entry.xgrid,
                        mode=mode,
                        xspace=xspace,
                        flavor_signs=perturb_signs,
                        calibration_gv=gv_flav_calib,
                    )
                    gv_pert_list.append(gv_pert)

                if sr_norm is not None:
                    gv_evol_list = obs.rotate_to_evolution(gv_pert_list)
                    gv_evol_list = [
                        analyzer._apply_norm_to_gv(
                            gv,
                            sr_norm,
                            range(14) if entry.hadronic else entry.flavor_indices,
                        )
                        for gv, entry in zip(gv_evol_list, obs.fk_entries)
                    ]
                    chi2_arr = obs.chi2(gv_evol_list)
                else:
                    chi2_arr = obs.chi2_from_flavor(gv_pert_list)
            else:
                gv_pert_list = []
                for idx, entry in enumerate(obs.fk_entries):
                    if entry.hadronic:
                        gv = analyzer._get_gv_all14_for_entry(obs, idx)
                        gv_calib = analyzer._get_calibration_gv_all14_for_entry(obs, idx)
                        perturb_players = list(coalition)
                        perturb_idx = [analyzer.flavor_indices[p] for p in coalition]
                    else:
                        gv = analyzer._get_gv_for_entry(obs, idx)
                        gv_calib = analyzer._get_calibration_gv_for_entry(obs, idx)
                        perturb_idx, perturb_players = (
                            analyzer._local_flavor_indices_and_players_for_entry(
                                entry, coalition
                            )
                        )
                    perturb_signs = (
                        sign_matrix[:, perturb_players]
                        if sign_matrix is not None else None
                    )

                    gv_pert = apply_gaussian_perturbation(
                        gv,
                        perturb_idx,
                        mu,
                        sigma,
                        amplitude,
                        entry.xgrid,
                        mode=mode,
                        xspace=xspace,
                        flavor_signs=perturb_signs,
                        calibration_gv=gv_calib,
                    )
                    if sr_norm is not None:
                        fi = range(14) if entry.hadronic else entry.flavor_indices
                        gv_pert = analyzer._apply_norm_to_gv(gv_pert, sr_norm, fi)
                    gv_pert_list.append(gv_pert)
                chi2_arr = obs.chi2(gv_pert_list)

            mean_samples.append(float(np.mean(chi2_arr)))

        mean_chi2 = float(np.mean(mean_samples))
        rows.append(
            {
                "dataset": obs.name,
                "ndata": int(obs.ndata),
                "operation": str(obs.operation),
                "n_fk": int(obs.n_fk),
                "hadronic": bool(obs.hadronic),
                "mean_chi2": mean_chi2,
                "chi2_per_point": mean_chi2 / obs.ndata,
            }
        )
    return rows


def _build_stabilization_report(analyzer, raw_results, pert, stab_cfg):
    """Build coalition->dataset outlier report and exclusion list."""
    diag = raw_results.get("diagnostic") or {}
    sign_matrices = raw_results.get("_sign_matrices")
    outliers = list(diag.get("outlier_coalitions", []))
    n_selected = min(len(outliers), int(stab_cfg["max_outlier_coalitions"]))
    selected = outliers[:n_selected]

    baseline_rows = _dataset_mean_chi2_for_coalition(
        analyzer, [], pert, sign_matrices=sign_matrices
    )
    baseline_map = {r["dataset"]: r for r in baseline_rows}

    threshold = float(stab_cfg["dataset_delta_chi2_threshold"])
    coalition_reports = []
    dataset_acc = {}

    for oc in selected:
        coalition_idx = list(oc.get("coalition", []))
        coalition_labels = list(oc.get("coalition_labels", []))
        coalition_rows = _dataset_mean_chi2_for_coalition(
            analyzer, coalition_idx, pert, sign_matrices=sign_matrices
        )

        flagged = []
        for crow in coalition_rows:
            dname = crow["dataset"]
            brow = baseline_map[dname]
            delta_mean = float(crow["mean_chi2"] - brow["mean_chi2"])
            if abs(delta_mean) < threshold:
                continue

            row = {
                "dataset": dname,
                "ndata": int(crow["ndata"]),
                "operation": crow["operation"],
                "n_fk": int(crow["n_fk"]),
                "hadronic": bool(crow["hadronic"]),
                "baseline_mean_chi2": float(brow["mean_chi2"]),
                "coalition_mean_chi2": float(crow["mean_chi2"]),
                "delta_mean_chi2": delta_mean,
                "baseline_chi2_per_point": float(brow["chi2_per_point"]),
                "coalition_chi2_per_point": float(crow["chi2_per_point"]),
                "delta_chi2_per_point": float(
                    crow["chi2_per_point"] - brow["chi2_per_point"]
                ),
            }
            flagged.append(row)

            acc = dataset_acc.setdefault(
                dname,
                {
                    "dataset": dname,
                    "ndata": int(crow["ndata"]),
                    "operation": crow["operation"],
                    "n_fk": int(crow["n_fk"]),
                    "hadronic": bool(crow["hadronic"]),
                    "max_abs_delta_mean_chi2": 0.0,
                    "coalitions": [],
                },
            )
            acc["max_abs_delta_mean_chi2"] = float(
                max(acc["max_abs_delta_mean_chi2"], abs(delta_mean))
            )
            acc["coalitions"].append(
                {
                    "coalition": coalition_idx,
                    "coalition_labels": coalition_labels,
                    "delta_mean_chi2": delta_mean,
                    "coalition_chi2": float(oc.get("chi2", float("nan"))),
                }
            )

        flagged.sort(key=lambda x: abs(x["delta_mean_chi2"]), reverse=True)
        coalition_reports.append(
            {
                "coalition": coalition_idx,
                "coalition_labels": coalition_labels,
                "size": int(oc.get("size", len(coalition_idx))),
                "chi2": float(oc.get("chi2", float("nan"))),
                "z_score": float(oc.get("z_score", float("nan"))),
                "n_flagged_datasets": len(flagged),
                "flagged_datasets": flagged,
            }
        )

    flagged_datasets = sorted(
        dataset_acc.values(),
        key=lambda x: abs(x["max_abs_delta_mean_chi2"]),
        reverse=True,
    )
    excluded = (
        [d["dataset"] for d in flagged_datasets]
        if stab_cfg["action"] == "exclude_dataset"
        else []
    )

    return {
        "enabled": bool(stab_cfg["enabled"]),
        "action": stab_cfg["action"],
        "dataset_delta_chi2_threshold": threshold,
        "max_outlier_coalitions": int(stab_cfg["max_outlier_coalitions"]),
        "n_outlier_coalitions_detected": len(outliers),
        "n_outlier_coalitions_analyzed": n_selected,
        "coalitions_analyzed": coalition_reports,
        "n_flagged_datasets": len(flagged_datasets),
        "flagged_datasets": flagged_datasets,
        "excluded_datasets": excluded,
    }


def _save_stabilization_files(report, output_dir, basis):
    """Write stabilization report JSON + compact flagged-dataset CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"stabilization_report_{basis}.json"
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Saved: {json_path}")

    csv_path = output_dir / f"stabilization_flagged_datasets_{basis}.csv"
    with open(csv_path, "w") as f:
        f.write(
            "dataset,operation,n_fk,ndata,max_abs_delta_mean_chi2,"
            "n_flagged_coalitions,excluded\n"
        )
        excluded = set(report.get("excluded_datasets", []))
        for row in report.get("flagged_datasets", []):
            f.write(
                f"{row['dataset']},{row['operation']},{row['n_fk']},{row['ndata']},"
                f"{row['max_abs_delta_mean_chi2']:.8f},{len(row.get('coalitions', []))},"
                f"{1 if row['dataset'] in excluded else 0}\n"
            )
    print(f"Saved: {csv_path}")
    return json_path, csv_path


def run_analysis(cfg, output_dir, setup_context, n_jobs_override=None,
                 diagnostic=None, outlier_n_sigma=3.0):
    """Execute the full Shapley value pipeline from one experiment config.

    Parameters
    ----------
    diagnostic : bool or None
        When True, record chi2 for every coalition and write diagnostic files.
        None (default) reads the ``diagnostic`` key from *cfg*.
    outlier_n_sigma : float
        Z-score threshold for flagging extreme coalitions / marginals.
    """
    pdf = setup_context["pdf"]
    observables = setup_context["observables"]
    flavor_info = setup_context["flavor_info"]

    n_replicas = cfg.get("n_replicas", 100)
    bases = cfg.get("basis", ["evolution"])
    if isinstance(bases, str):
        bases = [bases]

    pert = cfg["perturbation"]
    mu = pert["mu"]
    sigma = pert["sigma"]
    amplitude = pert["amplitude"]
    mode = pert.get("mode", "additive")
    xspace = pert.get("xspace", "linear")
    # Per-replica SVs: compute phi_j^(k) for every replica and report uncertainty.
    # Resolved from perturbation block or top-level cfg key.
    per_replica = bool(
        pert.get("per_replica", cfg.get("per_replica", False))
    )
    deterministic_sign_symmetrized = bool(
        pert.get(
            "deterministic_sign_symmetrized",
            cfg.get("deterministic_sign_symmetrized", True),
        )
    )
    random_sign = bool(
        pert.get("random_sign", cfg.get("random_sign", False))
    )
    n_sign_samples = int(
        pert.get("n_sign_samples", cfg.get("n_sign_samples", 1))
    )
    if deterministic_sign_symmetrized and random_sign:
        print(
            "Note: ignoring random_sign=True; using deterministic calibrated up/down games."
        )
    if deterministic_sign_symmetrized and n_sign_samples != 1:
        print(
            "Note: ignoring n_sign_samples="
            f"{n_sign_samples}; deterministic dual game uses one fixed definition."
        )
    random_seed = pert.get("random_seed", cfg.get("random_seed"))
    explicit_member_mode = pert.get("member_mode", cfg.get("member_mode"))
    if explicit_member_mode is not None:
        member_mode = str(explicit_member_mode).strip().lower()
        if member_mode not in {"replicas", "central"}:
            raise ValueError(
                "member_mode must be 'replicas' or 'central'; "
                f"got '{member_mode}'."
            )
    else:
        member_mode = "replicas" if per_replica else "central"
    enforce_sumrules = cfg.get("enforce_sumrules", False)
    n_jobs = int(cfg.get("n_jobs", 1))
    if n_jobs_override is not None:
        n_jobs = int(n_jobs_override)
    # Resolve diagnostic flag: CLI/caller override > runcard key > False.
    if diagnostic is None:
        diagnostic = bool(cfg.get("diagnostic", False))
    diag_sigma = float(cfg.get("outlier_n_sigma", outlier_n_sigma))

    stabilization = _resolve_stabilization_cfg(cfg)
    if stabilization["enabled"]:
        print(
            "Stabilization     : ON  "
            f"(action={stabilization['action']}, "
            f"threshold={stabilization['dataset_delta_chi2_threshold']:.1f}, "
            f"max_outlier_coalitions={stabilization['max_outlier_coalitions']})"
        )
    else:
        print("Stabilization     : OFF")

    output_dir.mkdir(parents=True, exist_ok=True)
    all_results = {}

    for basis in bases:
        print(f"\n{'=' * 60}")
        print(f"  Basis: {basis}")
        print(f"{'=' * 60}\n")

        if basis == "flavor":
            fi = {
                "indices": [2, 3, 4, 5, 6, 7, 8, 9, 10],
                "pdg_codes": [-4, -3, -2, -1, 21, 1, 2, 3, 4],
                "names": ["cbar", "sbar", "ubar", "dbar", "g", "d", "u", "s", "c"],
                "n_flavors": 9,
            }
        else:
            fi = flavor_info
        
        # Debugging
        #with open("/tmp/perturbation_debug.log", "a") as f:
        #            f.write(f"basis and mode:{basis} and {mode}\n")

        analyzer = NNPDFShapleyAnalyzer(
            pdf, observables, fi,
            n_replicas=n_replicas,
            basis=basis,
            enforce_sumrules=enforce_sumrules,
            member_mode=member_mode,
        )

        t0 = time.time()
        # Stabilization requires full coalition diagnostics from the raw run.
        raw_diagnostic = diagnostic
        if stabilization["enabled"]:
            if diagnostic is False:
                print(
                    "  Note: forcing diagnostics ON for raw run because "
                    "stabilization is enabled."
                )
            raw_diagnostic = True

        results_raw = analyzer.exact_shap(
            mu=mu, sigma=sigma, amplitude=amplitude,
            mode=mode, basis=basis, xspace=xspace, plot=True, n_jobs=n_jobs,
            diagnostic=raw_diagnostic, outlier_n_sigma=diag_sigma,
            per_replica=per_replica, random_sign=random_sign,
            n_sign_samples=n_sign_samples,
            random_seed=random_seed,
            deterministic_sign_symmetrized=deterministic_sign_symmetrized,
        )
        stabilization_report = None
        stabilization_json = None
        excluded_datasets = []
        results_final = results_raw
        stable_rerun_performed = False
        observables_used = observables

        if stabilization["enabled"]:
            stabilization_report = _build_stabilization_report(
                analyzer, results_raw, pert, stabilization
            )
            stabilization_json, _ = _save_stabilization_files(
                stabilization_report, output_dir, basis
            )
            excluded_datasets = list(stabilization_report.get("excluded_datasets", []))

            if (
                stabilization["action"] == "exclude_dataset"
                and stabilization["rerun_stable"]
                and excluded_datasets
            ):
                kept = [obs for obs in observables if obs.name not in excluded_datasets]
                if not kept:
                    print(
                        "  Stabilization requested exclusion, but all datasets would "
                        "be removed. Keeping raw result."
                    )
                else:
                    print(
                        f"  Excluding {len(excluded_datasets)} dataset(s) and "
                        f"re-running stable Shapley on {len(kept)} dataset(s)."
                    )
                    stable_analyzer = NNPDFShapleyAnalyzer(
                        pdf,
                        kept,
                        fi,
                        n_replicas=n_replicas,
                        basis=basis,
                        enforce_sumrules=enforce_sumrules,
                        member_mode=member_mode,
                    )
                    results_final = stable_analyzer.exact_shap(
                        mu=mu, sigma=sigma, amplitude=amplitude,
                        mode=mode, basis=basis, xspace=xspace, plot=True, n_jobs=n_jobs,
                        diagnostic=diagnostic, outlier_n_sigma=diag_sigma,
                        per_replica=per_replica, random_sign=random_sign,
                        n_sign_samples=n_sign_samples,
                        random_seed=random_seed,
                        deterministic_sign_symmetrized=deterministic_sign_symmetrized,
                    )
                    stable_rerun_performed = True
                    observables_used = kept

        elapsed = time.time() - t0
        print(f"\nElapsed: {elapsed:.1f}s")

        # Save plots from the final result (raw or stable rerun).
        if results_final.get("fig_pdfs") is not None:
            pdf_fig_path = output_dir / f"pdfs_{basis}.png"
            results_final["fig_pdfs"].savefig(pdf_fig_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {pdf_fig_path}")
        if results_final.get("fig_bar") is not None:
            bar_fig_path = output_dir / f"shapley_bar_{basis}.png"
            results_final["fig_bar"].savefig(bar_fig_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {bar_fig_path}")

        # Always keep raw reference output when stabilization is enabled.
        labels_raw = results_raw["flavor_short"]
        sv_raw = results_raw["shapley_values"]
        if stabilization["enabled"]:
            csv_raw = output_dir / f"shapley_values_{basis}_raw.csv"
            _write_shapley_csv(csv_raw, labels_raw, sv_raw)
            print(f"Saved: {csv_raw}")

        labels = results_final["flavor_short"]
        sv = results_final["shapley_values"]
        sv_odd = results_final.get("shapley_values_odd")
        csv_path = output_dir / f"shapley_values_{basis}.csv"
        _write_shapley_csv(csv_path, labels, sv)
        print(f"Saved: {csv_path}")
        if sv_odd is not None:
            csv_odd_path = output_dir / f"shapley_values_odd_{basis}.csv"
            _write_shapley_csv(csv_odd_path, labels, sv_odd)
            print(f"Saved: {csv_odd_path}")

        # Sampling-axis outputs.
        sv_per_replica = results_final.get("shapley_values_per_replica")
        sv_odd_per_replica = results_final.get("shapley_values_odd_per_replica")
        sv_odd_signed_per_replica = results_final.get("shapley_values_odd_signed_per_replica")
        sv_per_sample = results_final.get("shapley_values_per_sign_sample")
        sv_std = results_final.get("shapley_std")
        sv_err = results_final.get("shapley_err")
        sv_odd_std = results_final.get("shapley_odd_std")
        sv_odd_err = results_final.get("shapley_odd_err")
        sv_sign_std = results_final.get("shapley_sign_std")
        sv_sign_err = results_final.get("shapley_sign_err")
        if per_replica and sv_per_replica is not None:
            per_rep_path = output_dir / f"shapley_values_per_replica_{basis}.csv"
            _write_shapley_per_replica_csv(per_rep_path, labels, sv_per_replica)
            print(f"Saved: {per_rep_path}")
        if per_replica and sv_odd_per_replica is not None:
            per_rep_odd_path = output_dir / f"shapley_values_odd_per_replica_{basis}.csv"
            _write_shapley_per_replica_csv(
                per_rep_odd_path, labels, sv_odd_per_replica
            )
            print(f"Saved: {per_rep_odd_path}")
        if per_replica and sv_odd_signed_per_replica is not None:
            per_rep_odd_signed_path = (
                output_dir / f"shapley_values_odd_signed_per_replica_{basis}.csv"
            )
            _write_shapley_per_replica_csv(
                per_rep_odd_signed_path, labels, sv_odd_signed_per_replica
            )
            print(f"Saved: {per_rep_odd_signed_path}")
        if sv_per_sample is not None:
            per_sample_path = output_dir / f"shapley_values_per_sign_sample_{basis}.csv"
            _write_shapley_per_sample_csv(per_sample_path, labels, sv_per_sample)
            print(f"Saved: {per_sample_path}")
        if per_replica and sv_std is not None:
            unc_path = output_dir / f"shapley_uncertainties_{basis}.csv"
            with open(unc_path, "w") as f:
                f.write("flavour,mean,std,err\n")
                for lbl, mean_v, std_v, err_v in zip(
                    labels,
                    sv,
                    sv_std,
                    sv_err,
                ):
                    f.write(
                        f"{lbl},{float(mean_v):.8f},{float(std_v):.8f},{float(err_v):.8f}\n"
                    )
            print(f"Saved: {unc_path}")
        # Odd metric is the sign-sensitivity S_abs = <|phi_odd|>; no odd variance file.
        if sv_sign_std is not None:
            sign_unc_path = output_dir / f"shapley_sign_uncertainties_{basis}.csv"
            with open(sign_unc_path, "w") as f:
                f.write("flavour,mean,std,err\n")
                for lbl, mean_v, std_v, err_v in zip(
                    labels,
                    sv,
                    sv_sign_std,
                    sv_sign_err,
                ):
                    f.write(
                        f"{lbl},{float(mean_v):.8f},{float(std_v):.8f},{float(err_v):.8f}\n"
                    )
            print(f"Saved: {sign_unc_path}")

        if stable_rerun_performed:
            csv_stable = output_dir / f"shapley_values_{basis}_stable.csv"
            _write_shapley_csv(csv_stable, labels, sv)
            print(f"Saved: {csv_stable}")

        # Diagnostic files.
        if stabilization["enabled"] and results_raw.get("diagnostic"):
            _save_diagnostic_files(results_raw["diagnostic"], output_dir, f"{basis}_raw")
        if results_final.get("diagnostic"):
            _save_diagnostic_files(results_final["diagnostic"], output_dir, basis)

        all_results[basis] = {
            "shapley_values": {l: float(v) for l, v in zip(labels, sv)},
            "phi_even": {l: float(v) for l, v in zip(labels, sv)},
            "phi_odd": (
                {l: float(v) for l, v in zip(labels, sv_odd)}
                if sv_odd is not None else None
            ),
            "phi_odd_signed": (
                {l: float(v) for l, v in zip(labels, results_final.get("shapley_values_odd_signed", []))}
                if results_final.get("shapley_values_odd_signed") is not None else None
            ),
            "sign_sensitivity_abs": (
                {l: float(v) for l, v in zip(labels, sv_odd)}
                if sv_odd is not None else None
            ),
            "shapley_values_per_replica": _serialize_shapley_per_replica(
                labels, sv_per_replica
            ),
            "phi_even_per_replica": _serialize_shapley_per_replica(
                labels, sv_per_replica
            ),
            "phi_odd_per_replica": _serialize_shapley_per_replica(
                labels, sv_odd_per_replica
            ),
            "phi_odd_signed_per_replica": _serialize_shapley_per_replica(
                labels, sv_odd_signed_per_replica
            ),
            "shapley_values_per_sign_sample": _serialize_shapley_per_sample(
                labels, sv_per_sample
            ),
            "shapley_std": (
                {l: float(v) for l, v in zip(labels, sv_std)}
                if sv_std is not None else None
            ),
            "std_even_rep": (
                {l: float(v) for l, v in zip(labels, sv_std)}
                if sv_std is not None else None
            ),
            "shapley_err": (
                {l: float(v) for l, v in zip(labels, sv_err)}
                if sv_err is not None else None
            ),
            "shapley_values_odd": (
                {l: float(v) for l, v in zip(labels, sv_odd)}
                if sv_odd is not None else None
            ),
            "shapley_odd_std": (
                {l: float(v) for l, v in zip(labels, sv_odd_std)}
                if sv_odd_std is not None else None
            ),
            "std_odd_rep": (
                {l: float(v) for l, v in zip(labels, sv_odd_std)}
                if sv_odd_std is not None else None
            ),
            "shapley_odd_err": (
                {l: float(v) for l, v in zip(labels, sv_odd_err)}
                if sv_odd_err is not None else None
            ),
            "shapley_sign_std": (
                {l: float(v) for l, v in zip(labels, sv_sign_std)}
                if sv_sign_std is not None else None
            ),
            "shapley_sign_err": (
                {l: float(v) for l, v in zip(labels, sv_sign_err)}
                if sv_sign_err is not None else None
            ),
            "per_replica": per_replica,
            "random_sign": random_sign,
            "member_mode": member_mode,
            "n_sign_samples": int(results_final.get("n_sign_samples", n_sign_samples)),
            "antithetic_sign": bool(
                results_final.get("random_sign") and results_final.get("n_sign_samples", 1) > 1
            ),
            "deterministic_sign_symmetrized": bool(
                results_final.get(
                    "deterministic_sign_symmetrized",
                    deterministic_sign_symmetrized,
                )
            ),
            "random_seed": random_seed,
            "baseline_chi2": float(results_final["baseline_chi2"]),
            "validation_checks": results_final.get("checks"),
            "coalitions_evaluated": results_final["coalitions_evaluated"],
            "theory_evaluations": int(
                results_final.get("theory_evaluations", results_final["coalitions_evaluated"])
            ),
            "elapsed_seconds": round(elapsed, 1),
            "n_jobs": int(results_final.get("n_jobs", n_jobs)),
            "n_datasets_used": len(observables_used),
            "stabilization": {
                "enabled": bool(stabilization["enabled"]),
                "action": stabilization["action"],
                "dataset_delta_chi2_threshold": float(
                    stabilization["dataset_delta_chi2_threshold"]
                ),
                "max_outlier_coalitions": int(stabilization["max_outlier_coalitions"]),
                "stable_rerun_performed": bool(stable_rerun_performed),
                "n_excluded_datasets": int(len(excluded_datasets)),
                "excluded_datasets": excluded_datasets,
                "report_json": str(stabilization_json) if stabilization_json else None,
            },
        }

    meta = {
        "runcard": str(cfg.get("_source", "unknown")),
        "pdf_name": cfg["pdf_name"],
        "theory_id": cfg.get("theory_id", 708),
        "n_datasets": len(cfg["datasets"]),
        "n_replicas": n_replicas,
        "perturbation": {
            "mu": mu, "sigma": sigma, "amplitude": amplitude,
            "mode": mode, "xspace": xspace,
            "per_replica": per_replica,
            "random_sign": random_sign,
            "member_mode": member_mode,
            "n_sign_samples": int(n_sign_samples),
            "antithetic_sign": bool(random_sign and int(n_sign_samples) > 1),
            "deterministic_sign_symmetrized": bool(deterministic_sign_symmetrized),
            "random_seed": random_seed,
        },
        "enforce_sumrules": enforce_sumrules,
        "n_jobs": int(n_jobs),
        "stabilization": {
            "enabled": bool(stabilization["enabled"]),
            "action": stabilization["action"],
            "dataset_delta_chi2_threshold": float(
                stabilization["dataset_delta_chi2_threshold"]
            ),
            "max_outlier_coalitions": int(stabilization["max_outlier_coalitions"]),
            "rerun_stable": bool(stabilization["rerun_stable"]),
        },
    }
    combined = {"results_by_basis": all_results}
    json_path = str(output_dir / "results.json")
    save_results(combined, json_path, metadata=meta)
    print(f"\nSaved combined results: {json_path}")

    return all_results


def _save_consolidated_results(all_experiment_results, summary, output_dir):
    """Save consolidated results from all experiments into a single JSON.
    
    Combines results from all experiments into a single JSON file with a smart
    structure that preserves per-experiment metadata and results by basis.
    
    Parameters
    ----------
    all_experiment_results : dict
        Dictionary mapping experiment name to per-basis results
    summary : dict
        Experiment summary with perturbation and output_dir info
    output_dir : Path
        Top-level output directory
    """
    consolidated = {
        "timestamp": datetime.now().isoformat(),
        "experiments": {}
    }
    
    for exp_name, exp_results in all_experiment_results.items():
        exp_summary = summary.get(exp_name, {})
        consolidated["experiments"][exp_name] = {
            "output_dir": exp_summary.get("output_dir"),
            "perturbation": exp_summary.get("perturbation"),
            "results_by_basis": exp_results
        }
    
    json_path = output_dir / "consolidated_results.json"
    with open(json_path, "w") as f:
        json.dump(consolidated, f, indent=2)
    print(f"Saved consolidated results: {json_path}")


def main():
    parser = argparse.ArgumentParser(
        description="NNPDF Shapley value analysis from a YAML runcard."
    )
    parser.add_argument(
        "runcard", type=str,
        help="Path to YAML runcard."
    )
    parser.add_argument(
        "--output", "-o", type=str, default=None,
        help="Output directory. Default: sv_results/<runcard_stem>_<timestamp>."
    )
    parser.add_argument(
        "--n-jobs", type=int, default=None,
        help="Number of worker threads for coalition evaluation."
    )
    diag_group = parser.add_mutually_exclusive_group()
    diag_group.add_argument(
        "--diagnostic", action="store_true", default=None,
        help=(
            "Record chi2 for every coalition and write diagnostic files "
            "(coalition_chi2_<basis>.csv, marginal_contributions_<basis>.csv, "
            "diagnostic_stats_<basis>.json) into each experiment output directory. "
            "Overrides the 'diagnostic' key in the runcard."
        ),
    )
    diag_group.add_argument(
        "--no-diagnostic", dest="diagnostic", action="store_false",
        help="Disable diagnostics even if the runcard sets 'diagnostic: true'.",
    )
    parser.add_argument(
        "--outlier-n-sigma", type=float, default=3.0,
        help=(
            "Z-score threshold for flagging extreme coalitions and marginal "
            "contributions in the diagnostic output (default: 3.0)."
        ),
    )

    exp_sel = parser.add_argument_group("experiment selection")
    exp_sel.add_argument(
        "--list-experiments",
        action="store_true",
        help="List normalized experiments with indices and exit.",
    )
    exp_sel.add_argument(
        "--experiment",
        type=str,
        default=None,
        help=(
            "Run only the experiment with this name (as in the runcard). "
            "By default, single-experiment runs skip cross-experiment plots."
        ),
    )
    exp_sel.add_argument(
        "--experiment-index",
        type=int,
        default=None,
        help=(
            "Run only the experiment at this 0-based index in the normalized list. "
            "Useful for SLURM array tasks (SLURM_ARRAY_TASK_ID)."
        ),
    )
    exp_sel.add_argument(
        "--aggregate-only",
        action="store_true",
        help=(
            "Do not run Shapley; instead rebuild cross-experiment plots and "
            "consolidated summaries from per-experiment results already on disk."
        ),
    )
    exp_sel.add_argument(
        "--no-compare",
        action="store_true",
        help="Skip cross-experiment comparison plots and consolidated summaries.",
    )
    args = parser.parse_args()

    cfg = load_runcard(args.runcard)
    cfg["_source"] = args.runcard
    cfg = compute_ceilings(cfg)

    output_dir = _resolve_output_dir(cfg, args.runcard, args.output)

    print(f"Runcard : {args.runcard}")
    print(f"Output  : {output_dir}\n")

    experiments = _normalize_experiments(cfg)

    if args.list_experiments:
        print("Normalized experiments:")
        for i, (exp_name, exp_cfg) in enumerate(experiments):
            exp_token = _safe_token(exp_name) or "exp"
            pert = exp_cfg.get("perturbation", {})
            mu = pert.get("mu")
            sigma = pert.get("sigma")
            print(
                f"  [{i}] name={exp_name} token={exp_token} "
                f"mu={mu} sigma={sigma}"
            )
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve experiment selection (if any).
    selected = experiments
    selected_name = None
    if args.experiment_index is not None and args.experiment is not None:
        raise ValueError("Use only one of --experiment or --experiment-index")
    if args.experiment_index is not None:
        idx = int(args.experiment_index)
        if idx < 0 or idx >= len(experiments):
            raise IndexError(
                f"--experiment-index out of range: {idx} (n={len(experiments)})"
            )
        selected = [experiments[idx]]
        selected_name = experiments[idx][0]
    if args.experiment is not None:
        target = str(args.experiment)
        matches = [(n, c) for (n, c) in experiments if n == target]
        if not matches:
            # Convenience: allow matching by safe token as well.
            matches = [(n, c) for (n, c) in experiments if _safe_token(n) == _safe_token(target)]
        if not matches:
            known = ", ".join(n for (n, _) in experiments)
            raise KeyError(f"Unknown experiment '{target}'. Known: {known}")
        selected = [matches[0]]
        selected_name = matches[0][0]

    # By default, single-experiment runs skip cross-experiment comparison outputs.
    compare_enabled = not args.no_compare
    if len(selected) == 1 and not args.aggregate_only:
        compare_enabled = False

    def _load_results_json(path: Path):
        with open(path) as f:
            payload = json.load(f)
        results = payload.get("results") or {}
        rb = results.get("results_by_basis")
        if rb is None:
            raise KeyError(f"Missing results.results_by_basis in {path}")
        return rb, payload.get("metadata") or {}

    if args.aggregate_only:
        # Aggregation needs a deterministic output directory.
        if args.output is None and not cfg.get("output_dir"):
            raise ValueError(
                "--aggregate-only requires --output (or a runcard 'output_dir') "
                "so the timestamped run directory can be located."
            )

        summary = {}
        all_experiment_results = {}
        missing = []
        for exp_name, exp_cfg in experiments:
            exp_token = _safe_token(exp_name) or "exp"
            exp_output_dir = output_dir / exp_token
            res_path = exp_output_dir / "results.json"
            if not res_path.exists():
                missing.append(str(res_path))
                continue
            exp_results_by_basis, _ = _load_results_json(res_path)
            all_experiment_results[exp_name] = exp_results_by_basis
            summary[exp_name] = {
                "output_dir": str(exp_output_dir),
                "perturbation": exp_cfg["perturbation"],
            }

        if missing:
            print("WARNING: some experiments are missing results.json:")
            for p in missing:
                print(f"  - {p}")

        if compare_enabled:
            _save_comparison_plots(all_experiment_results, output_dir, summary)

            summary_path = output_dir / "experiments_summary.json"
            with open(summary_path, "w") as f:
                json.dump(summary, f, indent=2)
            print(f"\nSaved multi-experiment summary: {summary_path}")

            _save_consolidated_results(all_experiment_results, summary, output_dir)
        return

    setup_context = _build_setup_context(cfg)

    summary = {}
    all_experiment_results = {}
    for exp_name, exp_cfg in selected:
        exp_token = _safe_token(exp_name) or "exp"
        exp_output_dir = output_dir / exp_token
        #  For physical perturbation mode
        pert = exp_cfg['perturbation']
        if pert['mode'] == 'physical':
            perturbation.ceilings_plus = pert['ceilings_plus']
            perturbation.ceilings_minus = pert['ceilings_minus']
            print(f"\nCeilings updated for experiment '{exp_name}'")
        print(f"\nRunning experiment '{exp_name}'")
        print(f"Output  : {exp_output_dir}\n")
        experiment_results = run_analysis(
            exp_cfg,
            exp_output_dir,
            setup_context,
            n_jobs_override=args.n_jobs,
            diagnostic=args.diagnostic,
            outlier_n_sigma=args.outlier_n_sigma,
        )
        all_experiment_results[exp_name] = experiment_results
        summary[exp_name] = {
            "output_dir": str(exp_output_dir),
            "perturbation": exp_cfg["perturbation"],
        }

    if compare_enabled and len(all_experiment_results) > 0:
        _save_comparison_plots(all_experiment_results, output_dir, summary)

        summary_path = output_dir / "experiments_summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"\nSaved multi-experiment summary: {summary_path}")

        _save_consolidated_results(all_experiment_results, summary, output_dir)


if __name__ == "__main__":
    main()
