#!/usr/bin/env python3
"""Run NNPDF Shapley value analysis from a YAML runcard.

Usage
-----
    vp-shapley runcards/sv_dis_hera.yaml
    vp-shapley runcards/sv_global.yaml --output results/sv
"""

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml

from validphys.shapley.setup import setup_observables
from validphys.shapley.analyzer import NNPDFShapleyAnalyzer
from shapley_values import save_results, plot_shapley_comparison


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
        for key in ["basis", "n_jobs", "enforce_sumrules", "n_replicas", "run_name"]:
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


def _format_experiment_title(exp_name, exp_meta):
    """Build concise comparison subplot title including perturbation params."""
    pert = (exp_meta or {}).get("perturbation", {})
    amp = pert.get("amplitude", "na")
    mu = pert.get("mu", "na")
    sigma = pert.get("sigma", "na")
    mode = pert.get("mode", "additive")
    xspace = pert.get("xspace")
    if mode == "ablation":
        title = f"{exp_name} | {mode}"
    elif mode == "calibrated":
        title = f"{exp_name} | A={amp}\u03c3_rep, mu={mu}, sigma={sigma}, {mode}"
    else:
        title = f"{exp_name} | A={amp}, mu={mu}, sigma={sigma}, {mode}"
    if xspace is not None:
        title += f", {xspace}"
    return title


def _save_comparison_plots(all_experiment_results, output_dir, experiment_meta):
    """Save per-basis comparison plots across experiments."""
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

    for basis in basis_names:
        labels = None
        for exp_name in exp_names:
            if basis in all_experiment_results[exp_name]:
                labels = list(all_experiment_results[exp_name][basis]["shapley_values"].keys())
                break
        if labels is None:
            continue

        matrix = []
        for exp_name in exp_names:
            basis_results = all_experiment_results[exp_name].get(basis, {})
            sv_map = basis_results.get("shapley_values", {})
            matrix.append([float(sv_map.get(lbl, 0.0)) for lbl in labels])
        matrix = np.asarray(matrix, dtype=float)

        results_list = []
        for i, exp_name in enumerate(exp_names):
            results_list.append(
                {
                    "shapley_values": matrix[i],
                    "player_short": labels,
                }
            )
        titles = [
            _format_experiment_title(exp_name, experiment_meta.get(exp_name, {}))
            for exp_name in exp_names
        ]
        fig = plot_shapley_comparison(results_list, titles=titles)

        out_path = output_dir / f"shapley_comparison_{basis}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved comparison plot: {out_path}")


def run_analysis(cfg, output_dir, setup_context, n_jobs_override=None):
    """Execute the full Shapley value pipeline from one experiment config."""
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
    enforce_sumrules = cfg.get("enforce_sumrules", False)
    n_jobs = int(cfg.get("n_jobs", 1))
    if n_jobs_override is not None:
        n_jobs = int(n_jobs_override)

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

        analyzer = NNPDFShapleyAnalyzer(
            pdf, observables, fi,
            n_replicas=n_replicas,
            basis=basis,
            enforce_sumrules=enforce_sumrules,
        )

        t0 = time.time()
        results = analyzer.exact_shap(
            mu=mu, sigma=sigma, amplitude=amplitude,
            mode=mode, xspace=xspace, plot=True, n_jobs=n_jobs,
        )
        elapsed = time.time() - t0
        print(f"\nElapsed: {elapsed:.1f}s")

        if results.get("fig_pdfs") is not None:
            pdf_fig_path = output_dir / f"pdfs_{basis}.png"
            results["fig_pdfs"].savefig(pdf_fig_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {pdf_fig_path}")
        if results.get("fig_bar") is not None:
            bar_fig_path = output_dir / f"shapley_bar_{basis}.png"
            results["fig_bar"].savefig(bar_fig_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {bar_fig_path}")

        sv = results["shapley_values"]
        labels = results["flavor_short"]
        csv_path = output_dir / f"shapley_values_{basis}.csv"
        with open(csv_path, "w") as f:
            f.write("flavour,shapley_value\n")
            for lbl, val in zip(labels, sv):
                f.write(f"{lbl},{val:.8f}\n")
        print(f"Saved: {csv_path}")

        all_results[basis] = {
            "shapley_values": {l: float(v) for l, v in zip(labels, sv)},
            "baseline_chi2": float(results["baseline_chi2"]),
            "coalitions_evaluated": results["coalitions_evaluated"],
            "elapsed_seconds": round(elapsed, 1),
            "n_jobs": int(results.get("n_jobs", n_jobs)),
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
        },
        "enforce_sumrules": enforce_sumrules,
        "n_jobs": int(n_jobs),
    }
    combined = {"results_by_basis": all_results}
    json_path = str(output_dir / "results.json")
    save_results(combined, json_path, metadata=meta)
    print(f"\nSaved combined results: {json_path}")

    return all_results


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
    args = parser.parse_args()

    cfg = load_runcard(args.runcard)
    cfg["_source"] = args.runcard

    output_dir = _resolve_output_dir(cfg, args.runcard, args.output)

    print(f"Runcard : {args.runcard}")
    print(f"Output  : {output_dir}\n")

    experiments = _normalize_experiments(cfg)
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_context = _build_setup_context(cfg)

    summary = {}
    all_experiment_results = {}
    for exp_name, exp_cfg in experiments:
        exp_token = _safe_token(exp_name) or "exp"
        exp_output_dir = output_dir / exp_token
        print(f"\nRunning experiment '{exp_name}'")
        print(f"Output  : {exp_output_dir}\n")
        experiment_results = run_analysis(
            exp_cfg,
            exp_output_dir,
            setup_context,
            n_jobs_override=args.n_jobs,
        )
        all_experiment_results[exp_name] = experiment_results
        summary[exp_name] = {
            "output_dir": str(exp_output_dir),
            "perturbation": exp_cfg["perturbation"],
        }

    _save_comparison_plots(all_experiment_results, output_dir, summary)

    summary_path = output_dir / "experiments_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved multi-experiment summary: {summary_path}")


if __name__ == "__main__":
    main()
