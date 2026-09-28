"""
Resample pseudo-replicas from an already-trained BNN without retraining.
(DOOB-ed [Documented-Optimized-Organized-Beautified] by Claude Code (Sonnet 5))

Pipeline
--------
    vp-setupfit → n3fit (save=weights.weights.h5) → bnn_inference → evolven3fit → postfit

Use case: you trained N BNN models and drew M samples from each.  Now you want
M' > M samples — just reload the variational posterior and draw again.

Requirements
------------
The fit directory must contain, for each BNN model bnn_idx:

    nnfit/replica_{bnn_idx * n_bnn_samples_orig}/weights.weights.h5

i.e. n3fit must have been run with ``save: weights.weights.h5`` in the runcard
(or equivalently via the ``--save`` CLI flag).  Any pseudo-replica from the
same BNN model will do; the script tries all replicas in that model's block.

Both BNN sampling modes are supported post-hoc:
  - ``sampler="weight"``: fully self-contained, only needs the saved posterior
    weights (mu_w, logsig2_w, bias, flow weights).
  - ``sampler="function"`` (linearized Laplace): additionally needs the
    per-dataset training data (invcovmat, x-grid, FK-convolution
    ``ObservableWrapper``) that ``BNNPredictor.compute_sigma_theta`` uses to
    differentiate training predictions w.r.t. the MAP weights. That data is
    NOT embedded in the saved weights file, so it is reconstructed here by
    replaying the same reportengine/validphys data-loading resources
    (``n3fit_exec.py``'s own providers) that a normal training run resolves,
    stopping short of actually training the NN (see
    ``_reconstruct_training_data``). This needs the same theory/commondata
    to be available on disk that the original fit used (already true if you
    still have the runcard and haven't wiped your validphys/LHAPDF caches);
    it does not need the original fit's ``vp-setupfit`` output directory.

Usage
-----
    # Resample all BNN models in the fit with 100 samples each
    python bnn_inference.py sample --runcard path/to/runcard.yml --samples 100

    # Resample only SLURM task 3 (BNN model index 2) with 100 samples
    python bnn_inference.py sample --runcard path/to/runcard.yml --samples 100 --bnn-replica 3

    # Point at an explicit fit directory
    python bnn_inference.py sample --runcard runcard.yml --fit-dir /path/to/MyFit --samples 50

    # Force linearized-Laplace (function-space) resampling
    python bnn_inference.py sample --runcard runcard.yml --samples 100 --sampler function
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

log = logging.getLogger(__name__)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))


# ---------------------------------------------------------------------------
# Runcard parsing
# ---------------------------------------------------------------------------

def _load_runcard(runcard_path: Path) -> dict:
    """Return architecture and BNN metadata from an n3fit runcard."""
    import yaml

    with open(runcard_path) as fh:
        rc = yaml.safe_load(fh)

    params = rc["parameters"]
    return dict(
        nodes=params["nodes_per_layer"],
        activations=params["activation_per_layer"],
        initializer=params["initializer"],
        architecture=params["layer_type"],
        dropout_rate=params.get("dropout", 0.0),
        prior_prec=params.get("prior_prec", 0.1),
        std_init=params.get("std_init", None),
        dropout_rate_bayesian=params.get("bayes_dropout", 0.0),
        bayesian_bias=params.get("bayesian_bias", False),
        bayesian_flow=params.get("bayesian_flow", False),
        flow_n_couplings=params.get("flow_n_couplings", 3),
        flow_hidden_units=params.get("flow_hidden_units", 24),
        flav_info=rc["fitting"]["basis"],
        fitbasis=rc["fitting"]["fitbasis"],
        theoryid=rc["theory"]["theoryid"],
        n_bnn_models=params.get("n_bnn_models", 1),
        n_bnn_samples=params.get("n_bnn_samples", 1),
        sampling_space=params.get("sampling_space", "weight"),
    )


def _load_raw_runcard(runcard_path: Path) -> dict:
    """Return the full, unflattened runcard dict (needed by
    ``_reconstruct_training_data``, which resolves reportengine resources that
    expect the same nested structure -- ``dataset_inputs``, ``datacuts``,
    ``theory``, ``positivity``, ``integrability``, ``fitting`` -- n3fit itself
    parses via ``N3FitConfig``)."""
    import yaml

    with open(runcard_path) as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def _build_bnn_model(arch: dict, seed: int = 0):
    """Reconstruct a single-replica BNN pdf_model from architecture parameters."""
    from n3fit.model_gen import generate_pdf_model, ReplicaSettings

    replica_settings = ReplicaSettings(
        seed=seed,
        nodes=arch["nodes"],
        activations=arch["activations"],
        architecture=arch["architecture"],
        initializer=arch["initializer"],
        dropout_rate=arch["dropout_rate"],
        prior_prec=arch["prior_prec"],
        std_init=arch["std_init"],
        dropout_rate_bayesian=arch["dropout_rate_bayesian"],
        bayesian_bias=arch["bayesian_bias"],
        bayesian_flow=arch["bayesian_flow"],
        flow_n_couplings=arch["flow_n_couplings"],
        flow_hidden_units=arch["flow_hidden_units"],
    )
    return generate_pdf_model(
        replicas_settings=[replica_settings],
        flav_info=arch["flav_info"],
        fitbasis=arch["fitbasis"],
        impose_sumrule="All",
    )


def _reconstruct_training_data(rc_raw: dict, replica_cli: int):
    """
    Reconstruct the per-dataset training data (``invcovmat_per_dataset``,
    ``xgrid_per_dataset``, ``obs_wrappers_per_dataset``) that the linearized-Laplace
    sampler needs, WITHOUT retraining the NN.

    This works by resolving the same reportengine/validphys resources that a live
    n3fit run resolves as arguments to the ``performfit`` action --
    ``experiments_data``, ``replicas_nnseed_fitting_data_dict``,
    ``posdatasets_fitting_pos_dict``, ``integdatasets_fitting_integ_dict``, ``basis``,
    ``fitbasis`` (all defined in ``validphys.n3fit_data``/``validphys.results``,
    exposed via the same ``n3fit.scripts.n3fit_exec.N3FIT_PROVIDERS`` list) -- then
    feeding them into a real ``ModelTrainer`` and calling its own
    ``_generate_static_data``/``_generate_replica_losses`` (the exact methods
    ``hyperparametrizable`` uses before it starts training) to build the
    ``ObservableWrapper``s and invcovmats. No neural network is trained; only the
    (comparatively cheap, already-cached-on-disk) data-loading step runs.

    ``replica_cli`` must match the CLI replica number originally used to train this
    BNN model (``bnn_idx + 1``), since that is what seeds the train/validation mask
    generated for it.

    Parameters
    ----------
    rc_raw : dict
        The full, unflattened runcard dict (see ``_load_raw_runcard``).
    replica_cli : int
        1-indexed replica/SLURM-task number, matching what was passed to ``n3fit``
        on the command line when this BNN model was originally trained.

    Returns
    -------
    invcovmat_per_dataset, xgrid_per_dataset, obs_wrappers_per_dataset : list
        Per real-experimental-dataset (not positivity/integrability) data, in the
        same order, ready to be attached onto a rebuilt ``pdf_model``.
    """
    import importlib
    import pathlib

    from reportengine import namespaces, templateparser
    from reportengine.namespaces import NSList
    from reportengine.resourcebuilder import ResourceBuilder

    from n3fit.model_trainer import ModelTrainer
    from n3fit.scripts.n3fit_exec import N3FIT_PROVIDERS, N3FitConfig, N3FitEnvironment

    rc = dict(rc_raw)
    # Mirror N3FitConfig.from_yaml's fixed keys (n3fit_exec.py), on our own dict --
    # not the shared module-level N3FIT_FIXED_CONFIG, whose actions_ list is mutated
    # in place and must not accumulate across repeated calls/BNN models.
    rc["use_cuts"] = "internal"
    rc["use_t0"] = True
    rc["actions_"] = []
    rc["allow_legacy_names"] = False
    rc.setdefault("fiatlux", rc_raw.get("fiatlux"))
    rc.setdefault("positivity_bound", rc_raw.get("positivity_bound"))
    rc["use_thcovmat_in_fitting"] = False
    rc["use_thcovmat_in_sampling"] = False

    env = N3FitEnvironment()
    env.replicas = NSList([replica_cli], nskey="replica")
    env.hyperopt = None
    # Never written to: we only resolve data resources, we don't run the
    # `performfit` action (which is what would normally call init_output()).
    env.output_path = pathlib.Path(".") / f"_bnn_inference_scratch_{replica_cli}"
    env.replica_path = env.output_path / "nnfit"

    config = N3FitConfig(rc, environment=env)

    namespace = "datacuts::theory::fitting "
    resource_names = [
        "experiments_data",
        "replicas_nnseed_fitting_data_dict",
        "posdatasets_fitting_pos_dict",
        "integdatasets_fitting_integ_dict",
        "basis",
        "fitbasis",
    ]
    targets = [templateparser.string_to_target(namespace + name) for name in resource_names]
    providers = [importlib.import_module(p) for p in N3FIT_PROVIDERS]

    builder = ResourceBuilder(config, providers, targets, perform_final=False)
    builder.rootns.update(env.ns_dump())
    builder.resolve_fuzzytargets()
    builder.execute_sequential()
    ns = namespaces.resolve(builder.rootns, ("datacuts", "theory", "fitting"))

    experiments_data = ns["experiments_data"]
    replica, exp_info_single, nnseed = ns["replicas_nnseed_fitting_data_dict"][0]
    exp_info = [exp_info_single]
    pos_info = ns["posdatasets_fitting_pos_dict"]
    integ_info = ns["integdatasets_fitting_integ_dict"]
    basis = ns["basis"]
    fitbasis = ns["fitbasis"]

    model_trainer = ModelTrainer(
        experiments_data,
        exp_info,
        pos_info,
        integ_info,
        basis,
        fitbasis,
        [nnseed],
        rc.get("positivity_bound"),
        replicas=[replica],
    )
    # Only the observable/loss-generation half of hyperparametrizable -- no NN is
    # ever constructed or trained here.
    model_trainer._generate_static_data(None, None, None, None, 1, None)
    model_trainer._generate_replica_losses()

    return (
        model_trainer._experiment_data["invcovmat"],
        model_trainer._xgrid_per_dataset,
        model_trainer._obs_wrappers_per_dataset,
    )


# ---------------------------------------------------------------------------
# Weights discovery
# ---------------------------------------------------------------------------

def _detect_base_replica(nnfit_path: Path) -> int:
    """
    Return the lowest replica index under ``nnfit_path`` that has a saved
    ``weights.weights.h5``, i.e. BNN model 0's first pseudo-replica.

    Training always writes model ``bnn_idx``'s block as ``n_bnn_samples``
    *consecutive* replica directories, but the block for model 0 can start at
    replica 0 or replica 1 depending on how the fit was produced -- rather than
    assuming one or the other, detect it empirically so ``_find_source_weights``
    stays correct either way.
    """
    import re

    indices = [
        int(m.group(1))
        for p in nnfit_path.glob("replica_*/weights.weights.h5")
        if (m := re.match(r"replica_(\d+)$", p.parent.name))
    ]
    return min(indices) if indices else 0


def _find_source_weights(
    nnfit_path: Path, bnn_idx: int, n_bnn_samples_orig: int, base_replica: int = 0
):
    """
    Locate a saved weights file for BNN model ``bnn_idx``.

    During training, BNN model bnn_idx writes pseudo-replicas
    [base_replica + bnn_idx * n_bnn_samples_orig, ..., base_replica + (bnn_idx+1)
    * n_bnn_samples_orig - 1] (see ``_detect_base_replica``). All share the same
    variational posterior, so any of their weight files can be used to resample.
    """
    first = base_replica + bnn_idx * n_bnn_samples_orig
    for k in range(n_bnn_samples_orig):
        candidate = nnfit_path / f"replica_{first + k}" / "weights.weights.h5"
        if candidate.exists():
            return candidate
    return None


# ---------------------------------------------------------------------------
# Main sampling function
# ---------------------------------------------------------------------------

def sample_bnn(
    runcard,
    fit_dir=None,
    bnn_replica=None,
    samples: int = 100,
    sampler=None,
    start_replica=None,
):
    """
    Resample pseudo-replicas from a trained BNN posterior.

    Parameters
    ----------
    runcard : str or Path
        Path to the n3fit runcard YAML.
    fit_dir : str or Path, optional
        Fit output directory.  Defaults to a folder whose name matches
        the runcard stem, located next to the runcard.
    bnn_replica : int, optional
        1-indexed SLURM task ID of the BNN model to resample.
        If None, all n_bnn_models models are resampled in sequence.
    samples : int
        Number of pseudo-replicas to generate per BNN model.
    sampler : str, optional
        ``"weight"`` or ``"function"``.  Overrides the ``sampling_space``
        entry in the runcard.  ``"function"`` (linearized Laplace) additionally
        reconstructs the training data (see ``_reconstruct_training_data``);
        this needs the theory/commondata the original fit used to still be
        available on disk (validphys/LHAPDF caches), and is noticeably slower
        per BNN model than ``"weight"`` since it recomputes the Gauss-Newton
        posterior covariance.
    start_replica : int, optional
        First output replica index.  Defaults to ``bnn_idx * samples``
        so the new replicas occupy the same index block as the old ones.

    Returns
    -------
    int
        Total number of pseudo-replicas written.
    """
    import os
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

    runcard_path = Path(runcard).resolve()
    if not runcard_path.exists():
        raise FileNotFoundError(f"Runcard not found: {runcard_path}")

    arch = _load_runcard(runcard_path)
    rc_raw = _load_raw_runcard(runcard_path) if sampler == "function" or (
        sampler is None and arch["sampling_space"] == "function"
    ) else None
    sampler = sampler or arch["sampling_space"]

    n_bnn_models = arch["n_bnn_models"]
    n_bnn_samples_orig = arch["n_bnn_samples"]

    if fit_dir is None:
        fit_dir = runcard_path.parent / runcard_path.stem
    fit_dir = Path(fit_dir).resolve()
    nnfit_path = fit_dir / "nnfit"
    fitname = fit_dir.name

    if not nnfit_path.is_dir():
        raise FileNotFoundError(f"nnfit directory not found: {nnfit_path}")

    # Theory object needed by storefit for Q0 and quark thresholds
    from validphys.loader import Loader
    theory = Loader().check_theoryID(arch["theoryid"])
    q0 = theory.get_description().get("Q0")

    bnn_idx_list = [bnn_replica - 1] if bnn_replica is not None else list(range(n_bnn_models))

    from n3fit.bnn_wrapper_copy import BNNPredictor
    from n3fit.io.writer import storefit
    from n3fit.vpinterface import N3PDF

    # Detected once from disk rather than assumed, since BNN model 0's block can
    # start at replica 0 or replica 1 depending on how the fit was produced.
    base_replica = _detect_base_replica(nnfit_path)
    log.info(
        "Detected BNN model 0 starting at replica_%d (block size %d) under %s",
        base_replica, n_bnn_samples_orig, nnfit_path,
    )

    total_written = 0
    for bnn_idx in bnn_idx_list:
        weights_file = _find_source_weights(nnfit_path, bnn_idx, n_bnn_samples_orig, base_replica)
        if weights_file is None:
            first = base_replica + bnn_idx * n_bnn_samples_orig
            log.warning(
                "BNN model %d: no weights.weights.h5 found in replica_%d..%d. "
                "Ensure n3fit was run with save=weights.weights.h5 — skipping.",
                bnn_idx,
                first,
                first + n_bnn_samples_orig - 1,
            )
            continue

        log.info("BNN model %d: loading weights from %s", bnn_idx, weights_file)

        # Rebuild architecture and load variational posterior (mu_w, logsig2_w, bias,
        # and flow weights if bayesian_flow is set) -- this alone is enough for
        # sampler="weight".
        pdf_model = _build_bnn_model(arch, seed=bnn_idx)
        pdf_model.load_identical_replicas(str(weights_file))

        if sampler == "function":
            # sampler="function" additionally needs the per-dataset training data;
            # replica_cli = bnn_idx + 1 matches the CLI replica number originally
            # used to train this BNN model (see _reconstruct_training_data).
            log.info(
                "BNN model %d: reconstructing training data for linearized-Laplace "
                "sampling (this re-runs data loading, not training) ...",
                bnn_idx,
            )
            invcovmat, xgrid, obs_wrappers = _reconstruct_training_data(rc_raw, bnn_idx + 1)
            pdf_model.invcovmat_per_dataset = invcovmat
            pdf_model.xgrid_per_dataset = xgrid
            pdf_model.obs_wrappers_per_dataset = obs_wrappers

        # Sample from the posterior
        predictor = BNNPredictor(pdf_model=pdf_model, n_bnn_samples=samples, sampler=sampler)
        log.info("BNN model %d: generating %d %s-space samples ...", bnn_idx, samples, sampler)
        sampled_models = predictor.pdf_sampler()

        # Replica index block for output
        first_out = (start_replica + bnn_idx * samples) if start_replica is not None \
                    else bnn_idx * samples
        out_indices = list(range(first_out, first_out + samples))

        log.info(
            "BNN model %d: writing replicas %d..%d to %s",
            bnn_idx, out_indices[0], out_indices[-1], nnfit_path,
        )

        source_dir = weights_file.parent

        for pdf_model_s, rep_idx in zip(sampled_models, out_indices):
            out_dir = nnfit_path / f"replica_{rep_idx}"
            out_dir.mkdir(parents=True, exist_ok=True)

            n3pdf = N3PDF(pdf_model_s, fit_basis=arch["flav_info"], Q=q0)
            storefit(n3pdf, rep_idx, out_dir / f"{fitname}.exportgrid", theory)

            # Copy training metadata so postfit quality cuts have something to read.
            # These stats are from the BNN training run and identical across all
            # pseudo-replicas of the same model, which is correct.
            for fname in (f"{fitname}.json", "chi2exps.log"):
                src = source_dir / fname
                dst = out_dir / fname
                if src.exists() and not dst.exists():
                    shutil.copy2(src, dst)

            log.info("  replica_%d written", rep_idx)

        total_written += len(out_indices)
        log.info("BNN model %d: done.", bnn_idx)

    log.info(
        "Sampling complete. %d pseudo-replicas written to %s", total_written, nnfit_path
    )
    return total_written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description=(
            "Resample pseudo-replicas from a trained BNN without retraining.\n\n"
            "Requires the fit to have been run with save=weights.weights.h5 so "
            "that the variational posterior parameters are on disk."
        )
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    sp = subparsers.add_parser(
        "sample",
        help="Draw new pseudo-replicas from the saved BNN posterior",
    )
    sp.add_argument(
        "--runcard",
        type=str,
        required=True,
        help="Path to the n3fit runcard YAML",
    )
    sp.add_argument(
        "--samples",
        type=int,
        default=100,
        help="Pseudo-replicas to generate per BNN model (default: 100)",
    )
    sp.add_argument(
        "--fit-dir",
        type=str,
        default=None,
        help=(
            "Path to the fit output directory. "
            "Default: a folder named after the runcard stem, next to the runcard."
        ),
    )
    sp.add_argument(
        "--bnn-replica",
        type=int,
        default=None,
        help=(
            "SLURM array task ID (1-indexed) of the BNN model to resample. "
            "If omitted, all n_bnn_models models are resampled."
        ),
    )
    sp.add_argument(
        "--sampler",
        type=str,
        default=None,
        choices=["weight", "function"],
        help=(
            "Sampling method. Overrides the runcard's sampling_space. "
            "'weight': classical, self-contained. "
            "'function': linearized Laplace; reconstructs training data, slower."
        ),
    )
    sp.add_argument(
        "--start-replica",
        type=int,
        default=None,
        help=(
            "First output replica index. "
            "Default: bnn_idx * samples, matching the training-time numbering."
        ),
    )

    args = parser.parse_args()

    if args.mode == "sample":
        sample_bnn(
            runcard=args.runcard,
            fit_dir=args.fit_dir,
            bnn_replica=args.bnn_replica,
            samples=args.samples,
            sampler=args.sampler,
            start_replica=args.start_replica,
        )
