#!/usr/bin/env python3
"""Manually regenerate consolidated_results.json from existing sv_results."""

import json
from datetime import datetime
from pathlib import Path

def _save_consolidated_results(all_experiment_results, summary, output_dir):
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


def load_results_json(path: Path):
    with open(path) as f:
        payload = json.load(f)
    results = payload.get("results") or {}
    rb = results.get("results_by_basis")
    if rb is None:
        raise KeyError(f"Missing results.results_by_basis in {path}")
    return rb, payload.get("metadata") or {}


def main():
    sv_results = Path("sv_results")  # adjust if needed

    all_experiment_results = {}
    summary = {}
    missing = []

    for exp_dir in sorted(sv_results.iterdir()):
        if not exp_dir.is_dir():
            continue
        res_path = exp_dir / "results.json"
        if not res_path.exists():
            missing.append(str(res_path))
            continue

        results_by_basis, metadata = load_results_json(res_path)
        exp_name = exp_dir.name  # use folder name as experiment name

        all_experiment_results[exp_name] = results_by_basis
        summary[exp_name] = {
            "output_dir": str(exp_dir),
            "perturbation": metadata.get("perturbation"),
        }
        print(f"Loaded: {exp_dir.name}")

    if missing:
        print("\nWARNING: missing results.json in:")
        for p in missing:
            print(f"  - {p}")

    print(f"\nLoaded {len(all_experiment_results)} experiments. Generating consolidated JSON...")
    _save_consolidated_results(all_experiment_results, summary, sv_results)


if __name__ == "__main__":
    main()