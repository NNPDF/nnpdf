"""
MC Dropout inference for NNPDF dropout fits
======================================================================

Loads the trained weights of one replica, runs N stochastic forward passes
with dropout kept active (training=True), and computes the mean PDF over
those passes.  The result contains:

  - ``x``           : x-grid, shape (n_x,)
  - ``mean``        : mean PDF,  shape (n_x, 14)
  - ``std``         : std  PDF,  shape (n_x, 14)
  - ``samples``     : all N samples, shape (N, n_x, 14)
  - ``flavours``    : LHAPDF PID list, shape (14,)

Usage :
    python -m n3fit.mc_dropout_inference \\
        --fit-dir     nnpdf40-like-dropout-cluster \\
        --replica     1 \\
        --n-samples   200 \\
        --output      output/mc_dropout_central.npz

Requirements
------------
The fit directory must contain:
    nnfit/replica_<N>/weights.weights.h5   - saved Keras weights
    filter.yml or n3fit runcard            - architecture parameters

Warning: The architecture parameters are read from the n3fit runcard, not from the saved model.  
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

# ---------------------------------------------------------------------------
# Allow running as   python n3fit/src/n3fit/mc_dropout_inference.py   or as
#   python -m n3fit.mc_dropout_inference   (if src/ is on sys.path)
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))


# The pdf_model outputs 14 components in the n3fit EVOL evolution basis.
# These are NOT LHAPDF quark PIDs.  storefit() converts them to LHAPDF quark
# flavour order via evln2lha().  The .npz files produced here save the raw
# EVOL-basis output; mc_dropout_assemble.py applies evln2lha() when writing
# the .exportgrid files.
FLAVOUR_IDS = list(range(14))  # indices 0..13 for the 14 EVOL basis components

# Standard x-grid used by n3fit's writer.py (XGRID constant)
# We read it from writer.py to stay consistent; fall back to a built-in copy.
def _get_xgrid():
    try:
        from n3fit.io.writer import XGRID
        return XGRID
    except ImportError:
        pass
    # Fallback: identical grid used by n3fit/io/writer.py
    from n3fit.io.writer import XGRID  # noqa - re-raise if truly missing
    return XGRID

# Loading the architecture parameters from the fit runcard 
DEFAULT_RUNCARD = (
    Path(__file__).resolve().parent.parent.parent  # n3fit/src/../..  : n3fit
    / "runcards" / "examples" / "nnpdf40-like-dropout-cluster.yml"
)


def _load_architecture(runcard_path):
    """Return the architecture parameters needed to rebuild the pdf_model."""
    with open(runcard_path) as fh:
        rc = yaml.safe_load(fh)

    params = rc["parameters"]
    basis  = rc["fitting"]["basis"]
    fitbasis = rc["fitting"]["fitbasis"]

    return dict(
        nodes          = params["nodes_per_layer"],        
        activations    = params["activation_per_layer"],  
        initializer    = params["initializer"],            
        architecture   = params["layer_type"],             
        dropout_rate   = params.get("dropout", 0.0),      
        flav_info      = basis,
        fitbasis       = fitbasis,                         
    )


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def build_pdf_model(arch, seed=0):
    """
    Reconstruct a single-replica pdf_model with the given architecture.

    Parameters
    ----------
    arch : dict
        Output of ``_load_architecture``.
    seed : int
        Seed for the replica (affects weight initialisation, which is then
        overwritten when you call ``load_weights``).

    Returns
    -------
    pdf_model : MetaModel
    """
    from n3fit.model_gen import generate_pdf_model, ReplicaSettings

    replica_settings = ReplicaSettings(
        seed         = seed,
        nodes        = arch["nodes"],
        activations  = arch["activations"],
        architecture = arch["architecture"],
        initializer  = arch["initializer"],
        dropout_rate = arch["dropout_rate"],
    )

    pdf_model = generate_pdf_model(
        replicas_settings = [replica_settings],
        flav_info         = arch["flav_info"],
        fitbasis          = arch["fitbasis"],
        impose_sumrule    = "All",
    )
    return pdf_model


# ---------------------------------------------------------------------------
# Weights loading
# ---------------------------------------------------------------------------

def load_weights(pdf_model, weights_path):
    """
    Load Keras weights from ``weights_path`` into ``pdf_model``.

    The saved file is a single-replica model so we use
    ``load_identical_replicas``, which goes through the
    ``single_replica_generator`` to match layer names correctly.

    Parameters
    ----------
    pdf_model : MetaModel
    weights_path : str or Path
    """
    pdf_model.load_identical_replicas(str(weights_path))


# ---------------------------------------------------------------------------
# MC Dropout inference
# ---------------------------------------------------------------------------

def mc_dropout_central_value(pdf_model, x, n_samples=100):
    """
    Run ``n_samples`` stochastic forward passes with dropout active.

    Parameters
    ----------
    pdf_model : MetaModel
        Model with weights already loaded.
    x : np.ndarray
        x-grid, shape (n_x,).  Will be reshaped to (1, n_x, 1) for the model.
    n_samples : int
        Number of MC Dropout samples.

    Returns
    -------
    samples : np.ndarray, shape (n_samples, n_x, 14)
    mean    : np.ndarray, shape (n_x, 14)
    std     : np.ndarray, shape (n_x, 14)
    """
    # Model expects input shape (1, n_x, 1) : batch=1, n_x points, 1 feature (x)
    x_input = x.reshape(1, -1, 1).astype(np.float64)

    # mc_dropout_predict returns (all_samples, mean, std)
    # all_samples shape: (n_samples, 1, 1, n_x, 14)  [batch, n_replicas, n_x, 14]
    raw_samples, raw_mean, raw_std = pdf_model.mc_dropout_predict(
        x={"pdf_input": x_input},
        n_samples=n_samples,
    )

    # Squeeze batch (1) and replica (1) axes : (n_samples, n_x, 14)
    samples = raw_samples[:, 0, 0, :, :]   # (N, n_x, 14)
    mean    = raw_mean[0, 0, :, :]          # (n_x, 14)
    std     = raw_std[0, 0, :, :]           # (n_x, 14)

    return samples, mean, std


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_results(output_path, x, samples, mean, std):
    """Save results to a compressed numpy archive."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output_path,
        x        = x,
        samples  = samples,           # (N, n_x, 14)
        mean     = mean,              # (n_x, 14)
        std      = std,               # (n_x, 14)
        flavours = np.array(FLAVOUR_IDS),
    )
    print(f"Saved MC Dropout results : {output_path}")


def print_summary(x, mean, std, n_samples):
    """Print a brief summary of the inference results."""
    print(f"\n{'='*60}")
    print(f"MC Dropout Inference : Option A: Central Value")
    print(f"{'='*60}")
    print(f"  Number of stochastic samples : {n_samples}")
    print(f"  x-grid points                : {len(x)}")
    print(f"  x range                      : [{x.min():.2e}, {x.max():.2e}]")
    print(f"  Output shape                 : {mean.shape}  (n_x, 14 EVOL basis components)")
    print()
    print(f"  EVOL basis component indices: {FLAVOUR_IDS}")
    print()
    # Show gluon (index 6 in EVOL basis) and u-quark (index 8) at a few x values
    # Note: output is in EVOL basis, not LHAPDF quark flavours. mc_dropout_assemble.py
    # will apply evln2lha() transformation before writing .exportgrid files.
    gluon_idx = 6   # gluon in EVOL basis
    uquark_idx = 8  # u-quark in EVOL basis
    x_show    = [1e-4, 1e-2, 0.1, 0.3]
    print(f"  Sample PDF values (mean +/- std) in EVOL basis:")
    print(f"  {'x':>10}  {'component 6 (G)':>18}  {'component 8 (u)':>18}")
    for xs in x_show:
        idx = np.argmin(np.abs(x - xs))
        val6  = x[idx] * mean[idx, gluon_idx]
        std6  = x[idx] * std[idx, gluon_idx]
        val8  = x[idx] * mean[idx, uquark_idx]
        std8  = x[idx] * std[idx, uquark_idx]
        print(f"  {x[idx]:>10.3e}  {val6:>10.4f} +/- {std6:.4f}  {val8:>10.4f} +/- {std8:.4f}")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="MC Dropout inference (Option A: central value) for a dropout NNPDF fit"
    )
    parser.add_argument(
        "--fit-dir",
        default="nnpdf40-like-dropout-cluster",
        help="Path to the fit directory (contains nnfit/). Default: nnpdf40-like-dropout-cluster",
    )
    parser.add_argument(
        "--replica",
        type=int,
        default=1,
        help="Replica number to load weights from. Default: 1",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of MC Dropout forward passes. Default: 100",
    )
    parser.add_argument(
        "--runcard",
        default=None,
        help=(
            "Path to the n3fit runcard YAML. "
            f"Default: {DEFAULT_RUNCARD}"
        ),
    )
    parser.add_argument(
        "--output",
        default="output/mc_dropout_central.npz",
        help="Output numpy archive path. Default: output/mc_dropout_central.npz",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed for model initialisation (overwritten by loaded weights). Default: 0",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    fit_dir   = Path(args.fit_dir)
    weights_path = fit_dir / "nnfit" / f"replica_{args.replica}" / "weights.weights.h5"

    if not weights_path.exists():
        sys.exit(f"ERROR: weights file not found: {weights_path}")

    runcard_path = Path(args.runcard) if args.runcard else DEFAULT_RUNCARD
    if not runcard_path.exists():
        sys.exit(f"ERROR: runcard not found: {runcard_path}")

    print(f"\nRuncard   : {runcard_path}")
    print(f"Weights   : {weights_path}")
    print(f"MC samples: {args.n_samples}")

    # 1. Parse architecture from runcard
    print("\n[1/4] Parsing architecture from runcard ...")
    arch = _load_architecture(runcard_path)
    print(f"  nodes      = {arch['nodes']}")
    print(f"  activations= {arch['activations']}")
    print(f"  dropout    = {arch['dropout_rate']}")
    print(f"  fitbasis   = {arch['fitbasis']}")

    # 2. Build pdf_model
    print("\n[2/4] Building PDF model ...")
    # Suppress TF info messages during model build
    import os; os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    pdf_model = build_pdf_model(arch, seed=args.seed)
    print(f"  Model built: {pdf_model.name}")

    # 3. Load weights
    print(f"\n[3/4] Loading weights from replica {args.replica} ...")
    load_weights(pdf_model, weights_path)
    print("  Weights loaded.")

    # 4. Run MC Dropout
    print(f"\n[4/4] Running {args.n_samples} MC Dropout forward passes ...")
    from n3fit.io.writer import XGRID as x_grid
    samples, mean, std = mc_dropout_central_value(
        pdf_model, x_grid, n_samples=args.n_samples
    )

    # Summary + save
    print_summary(x_grid, mean, std, args.n_samples)
    save_results(args.output, x_grid, samples, mean, std)


if __name__ == "__main__":
    main()
