# vp-shapley

Compute Shapley values to quantify the relative importance of different PDF regions to global fit quality. This tool analyzes how perturbations in specific kinematic regions affect the chi-squared fit, providing interpretable metrics for feature importance in PDF fits.

## Installation

Requires the external `shapley-values` package:

    pip install git+https://github.com/rbonnetguerrini/shapley-values.git@main

## Quick Start

Local run (single machine):

    vp-shapley runcards/DIS/HERA.yaml --output results/my_analysis

Cluster run with parallelization:

    sbatch --cpus-per-task=64 --export=ALL,N_JOBS=64 \
      /path/to/slurm/run_vp_shapley.slurm \
      /path/to/runcards/DIS/HERA.yaml

## Usage

    vp-shapley <runcard.yaml> [OPTIONS]

    Options:
      --output, -o PATH        Output directory (default: sv_results/<name>_<timestamp>)
      --n-jobs N               Parallel worker threads for coalition evaluation
      --diagnostic             Write per-coalition chi2 and marginal contribution files
      --no-diagnostic          Disable diagnostics even if runcard sets diagnostic: true
      --outlier-n-sigma FLOAT  Z-score threshold for flagging extreme coalitions (default: 3.0)

## Runcard keys

    pdf_name         LHAPDF set name                      [required]
    datasets         List of dataset names                 [required]
    experiments      List / mapping of experiment configs  [required]
      name           Experiment label
      perturbation:
        mu           Perturbation centre in x
        sigma        Gaussian width
        amplitude    Perturbation amplitude
        mode         additive (default) | multiplicative
        xspace       linear (default) | logx
    theory_id        NNPDF theory ID                       [default: 708]
    basis            evolution (default) | flavor | list of both
    n_replicas       Number of PDF replicas                [default: 100]
    n_jobs           Worker threads (overridden by --n-jobs)
    enforce_sumrules true | false                          [default: false]
    per_replica      Report per-replica SV uncertainty     [default: false]
    random_sign      Flip amplitude sign per replica       [default: false]
    diagnostic       Save coalition-level diagnostics      [default: false]
    outlier_n_sigma  Z-score for outlier flagging          [default: 3.0]
    stabilization:
      enabled                     [default: true]
      action                      exclude_dataset | report_only
      dataset_delta_chi2_threshold[default: 1e4]
      max_outlier_coalitions      [default: 5]
      rerun_stable                Re-run after exclusion  [default: true]
    output_dir       Fixed output path (relative to shapley root or absolute)
    output_root      Root for timestamped output dirs
    run_name         Label used in output directory name

## Output (per experiment sub-directory)

**Main results:**
- `results.json` — Combined Shapley values and metadata
- `shapley_values_<basis>.csv` — Mean Shapley value per flavour
- `shapley_bar_<basis>.png` — Visualization of feature importance

**Uncertainties & Stability:**
- `shapley_uncertainties_<basis>.csv` — Per-replica uncertainties (if `per_replica: true`)
- `shapley_values_<basis>_raw.csv` — Pre-stabilization result
- `shapley_values_<basis>_stable.csv` — Post-exclusion result (if stabilization applied)
- `stabilization_report_<basis>.json` — Flagged datasets and exclusion reasoning

**Diagnostics (with `--diagnostic`):**
- `coalition_chi2_<basis>.csv` — Chi-squared per coalition
- `marginal_contributions_<basis>.csv` — Marginal value contributions
- `diagnostic_stats_<basis>.json` — Statistical summary

**Multi-experiment output** (when running multiple experiments):

    shapley_comparison_<basis>.png/pdf        SV vs scan parameter (binned)
    shapley_comparison_truex_<basis>.png/pdf  SV vs scan parameter (true values)
    shapley_comparison_bars_<basis>.png/pdf   Multi-panel bar comparison
    experiments_summary.json

## Interpreting Results

- **Positive SV**: Region is helpful to fit quality; removing it would increase global chi-squared
- **Negative SV**: Region conflicts with fit; perturbing it actually improves overall fit
- **Large |SV|**: Region is critical for PDF determination
- **Small SV**: Region has minimal impact on fit quality

## Stabilization

When enabled, the tool automatically detects and excludes datasets that cause numerical instability (e.g., extreme outlier coalitions). The `stabilization_report_<basis>.json` logs which datasets were flagged and why. Disable with `stabilization: {enabled: false}` or inspect flagged datasets before re-running.
