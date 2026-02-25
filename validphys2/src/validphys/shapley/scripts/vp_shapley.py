#!/usr/bin/env python3
"""Run NNPDF Shapley value analysis from a YAML runcard.

Usage
-----
    vp-shapley runcards/sv_dis_hera.yaml
    vp-shapley runcards/sv_global.yaml --output results/sv
"""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import yaml

from validphys.shapley.setup import setup_observables
from validphys.shapley.analyzer import NNPDFShapleyAnalyzer
from shapley_values import save_results


def load_runcard(path):
    """Load and validate a YAML runcard."""
    with open(path) as f:
        cfg = yaml.safe_load(f)

    required = ["pdf_name", "datasets", "perturbation"]
    for key in required:
        if key not in cfg:
            raise KeyError(f"Missing required key '{key}' in runcard")

    pert = cfg["perturbation"]
    for key in ["mu", "sigma", "amplitude"]:
        if key not in pert:
            raise KeyError(f"Missing perturbation.{key} in runcard")

    return cfg


def run_analysis(cfg, output_dir):
    """Execute the full Shapley value pipeline from a parsed runcard."""
    pdf, observables, flavor_info, flavor_basis_info = setup_observables(
        pdf_name=cfg["pdf_name"],
        datasets=cfg["datasets"],
        theory_id=cfg.get("theory_id", 708),
        use_cuts=cfg.get("use_cuts", "internal"),
        variant=cfg.get("variant", None),
    )

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
                "names": ["cbar", "sbar", "ubar", "dbar",
                          "g", "d", "u", "s", "c"],
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
            mode=mode, xspace=xspace, plot=True,
        )
        elapsed = time.time() - t0
        print(f"\nElapsed: {elapsed:.1f}s")

        # Save plots
        if results.get("fig_pdfs") is not None:
            pdf_fig_path = output_dir / f"pdfs_{basis}.png"
            results["fig_pdfs"].savefig(
                pdf_fig_path, dpi=150, bbox_inches="tight"
            )
            print(f"Saved: {pdf_fig_path}")
        if results.get("fig_bar") is not None:
            bar_fig_path = output_dir / f"shapley_bar_{basis}.png"
            results["fig_bar"].savefig(
                bar_fig_path, dpi=150, bbox_inches="tight"
            )
            print(f"Saved: {bar_fig_path}")

        # Save CSV
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
        }

    # Save combined JSON using shapley_values.save_results
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
    args = parser.parse_args()

    cfg = load_runcard(args.runcard)
    cfg["_source"] = args.runcard

    if args.output:
        output_dir = Path(args.output)
    else:
        stem = Path(args.runcard).stem
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path("sv_results") / f"{ts}_{stem}"

    print(f"Runcard : {args.runcard}")
    print(f"Output  : {output_dir}\n")

    run_analysis(cfg, output_dir)


if __name__ == "__main__":
    main()
