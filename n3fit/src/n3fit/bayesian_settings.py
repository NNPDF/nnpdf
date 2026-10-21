"""
bayesian_settings.py

All options of a Bayesian-NN fit live in one block of the runcard, ``parameters::bayesian``
(the architecture itself, ``parameters::layer_type`` with ``VBDense`` layers, stays where it is):

    parameters:
      layer_type: [dense, dense, VBDense]
      bayesian:                     # values shown = the defaults of a BNN fit (BAYESIAN_DEFAULTS)
        n_bnn_models: 1             # trained BNN models (one per SLURM task)
        n_bnn_samples: 1200         # posterior samples written per model
        sampling_space: weight      # weight | function (linearized Laplace)
        prior_prec: 1.0
        std_init: -7.0              # initial log-variance of the Bayesian weights
        bayesian_bias: [false, false, true]   # default: true on the Bayesian (VBDense) layers
        bayesian_preproc: false
        bayes_dropout: 0.0
        bayesian_flow: false        # plus flow_n_couplings: 3, flow_hidden_units: 24
        rank: 4                     # VBDense_correlated only, plus u_init: 1e-4
        kl_annealing: linear        # KL warm-up: linear (0 -> 1) | cosine (0.1 -> 1)
        kl_warmup_epochs: 3000      # length of the KL warm-up (default 3000, at most epochs // 4)
        train_sampling: coherent    # coherent: one weight draw per training step, shared by all x
                                    #   points and all calls of the layer | lrt: local
                                    #   reparameterization trick (independent noise per x point)
        per_replica:                # per-sample quantities instead of copies of the trained model's
          chi2: true                # chi2, erf_tr, erf_vl of each sample (postfit's chi2 veto then acts per sample)
          positivity: true          # pos_state of each sample, with n3fit's own stopping positivity check
          postfit_veto: true        # postfit: explicit per-replica positivity veto (recomputed from the
                                    #   exportgrid for replicas that were not checked by n3fit); also
                                    #   enables postfit's filling (chi2 <= 2x the passing mean)
          posterior_stopping_samples: 0   # > 0: stopping requires positivity at the posterior mean
                                          #   and at this many fixed posterior draws
          rejection_sampling: false       # keep only posterior draws that pass positivity
          rejection_max_draws: null       # default 50 * n_bnn_samples

Runcards that define these keys directly under ``parameters`` (the format used so far) keep
working: those values are read when the ``bayesian`` block does not define them. Keys set nowhere
take the defaults above, for BNN fits only (``layer_type`` with VBDense / VBDense_correlated).

Other posteriors (set these explicitly):
  - lrt training (as before 1ee4b2551): train_sampling: lrt, std_init: 0.0 (= -log(prior_prec), the
    prior width; LRT trains from there), kl_annealing: cosine, kl_warmup_epochs: epochs // 4 (the old
    repeating cosine warm-up no longer exists); not available with bayesian_flow.
  - normalizing flow on a Bayesian layer: bayesian_flow: [false, false, true] (per layer, like
    bayesian_bias), flow_n_couplings: 3, flow_hidden_units: 24; train_sampling: coherent (lrt is not
    implemented for flows), std_init: -7.0.
  - correlated low-rank posterior: layer_type [..., VBDense_correlated], rank: 4, u_init: 1e-4;
    train_sampling: coherent, std_init: -7.0.

This module has no backend imports so that it can be used by scripts that only parse runcards.
"""

# Bayesian options and the per-sample switches, with their defaults
BAYESIAN_KEYS = (
    "n_bnn_models",
    "n_bnn_samples",
    "sampling_space",
    "prior_prec",
    "std_init",
    "bayesian_bias",
    "bayesian_preproc",
    "bayes_dropout",
    "bayesian_flow",
    "flow_n_couplings",
    "flow_hidden_units",
    "rank",
    "u_init",
    "kl_annealing",
    "kl_warmup_epochs",
    "train_sampling",
)
KL_ANNEALING_SCHEDULES = ("linear", "cosine")
TRAIN_SAMPLING_MODES = ("coherent", "lrt")
BAYESIAN_LAYER_TYPES = ("VBDense", "VBDense_correlated")
# Defaults of a BNN fit for the keys that are set nowhere (bayesian_bias and kl_warmup_epochs depend
# on the rest of the runcard, see bayesian_parameters)
BAYESIAN_DEFAULTS = {
    "n_bnn_models": 1,
    "n_bnn_samples": 1200,
    "sampling_space": "weight",
    "prior_prec": 1.0,
    "std_init": -7.0,
    "bayesian_preproc": False,
    "bayes_dropout": 0.0,
    "bayesian_flow": False,
    "flow_n_couplings": 3,
    "flow_hidden_units": 24,
    "rank": 4,
    "u_init": 1e-4,
    "kl_annealing": "linear",
    "train_sampling": "coherent",
}
KL_WARMUP_DEFAULT = 3000
PER_REPLICA_DEFAULTS = {
    "chi2": True,
    "positivity": True,
    "postfit_veto": True,
    "posterior_stopping_samples": 0,
    "rejection_sampling": False,
    "rejection_max_draws": None,
}


class BayesianSettingsError(ValueError):
    pass


def check_bayesian_parameters(parameters):
    """Validate ``parameters::bayesian``: known keys only, no conflicting values with keys
    given directly under ``parameters``, and sensible per-replica switches."""
    block = parameters.get("bayesian")
    if block is None:
        return
    if not isinstance(block, dict):
        raise BayesianSettingsError("parameters::bayesian must be a mapping")
    unknown = set(block) - set(BAYESIAN_KEYS) - {"per_replica"}
    if unknown:
        raise BayesianSettingsError(
            f"Unknown key(s) in parameters::bayesian: {sorted(unknown)}. "
            f"Allowed: {list(BAYESIAN_KEYS) + ['per_replica']}"
        )
    for key in BAYESIAN_KEYS:
        if key in block and key in parameters and parameters[key] != block[key]:
            raise BayesianSettingsError(
                f"'{key}' is set both in parameters ({parameters[key]}) and in "
                f"parameters::bayesian ({block[key]}); keep only the one in parameters::bayesian"
            )
    if block.get("kl_annealing", "linear") not in KL_ANNEALING_SCHEDULES:
        raise BayesianSettingsError(
            f"parameters::bayesian::kl_annealing must be one of {KL_ANNEALING_SCHEDULES}"
        )
    warmup = block.get("kl_warmup_epochs")
    if warmup is not None and (isinstance(warmup, bool) or not isinstance(warmup, int) or warmup < 1):
        raise BayesianSettingsError("parameters::bayesian::kl_warmup_epochs must be a positive integer")
    if block.get("train_sampling", "coherent") not in TRAIN_SAMPLING_MODES:
        raise BayesianSettingsError(
            f"parameters::bayesian::train_sampling must be one of {TRAIN_SAMPLING_MODES}"
        )
    if block.get("train_sampling") == "lrt" and block.get("bayesian_flow"):
        raise BayesianSettingsError("train_sampling: lrt is not implemented for bayesian_flow")
    per_replica = block.get("per_replica") or {}
    if not isinstance(per_replica, dict):
        raise BayesianSettingsError("parameters::bayesian::per_replica must be a mapping")
    unknown = set(per_replica) - set(PER_REPLICA_DEFAULTS)
    if unknown:
        raise BayesianSettingsError(
            f"Unknown key(s) in parameters::bayesian::per_replica: {sorted(unknown)}. "
            f"Allowed: {list(PER_REPLICA_DEFAULTS)}"
        )
    for key in ("chi2", "positivity", "postfit_veto", "rejection_sampling"):
        if key in per_replica and not isinstance(per_replica[key], bool):
            raise BayesianSettingsError(f"per_replica::{key} must be true or false")
    samples = per_replica.get("posterior_stopping_samples", 0)
    if not isinstance(samples, int) or samples < 0:
        raise BayesianSettingsError("per_replica::posterior_stopping_samples must be an integer >= 0")
    max_draws = per_replica.get("rejection_max_draws")
    if max_draws is not None and (not isinstance(max_draws, int) or max_draws < 1):
        raise BayesianSettingsError("per_replica::rejection_max_draws must be a positive integer")


def is_bnn(parameters):
    """True if the architecture has a Bayesian layer"""
    # layer_type can arrive as a numpy array (standard fits), whose truth value is ambiguous
    layer_type = parameters.get("layer_type")
    if layer_type is None:
        return False
    return any(str(layer) in BAYESIAN_LAYER_TYPES for layer in layer_type)


def bayesian_parameters(parameters):
    """
    The Bayesian options of a fit, from ``parameters::bayesian`` with a fallback to the same
    keys directly under ``parameters``. For a BNN fit, keys set nowhere take the defaults
    (BAYESIAN_DEFAULTS; bayesian_bias: true on the Bayesian layers; kl_warmup_epochs:
    KL_WARMUP_DEFAULT, at most epochs // 4). For other fits they are absent, so that callers keep
    their own defaults. ``per_replica`` is always present, filled with ``PER_REPLICA_DEFAULTS``.
    """
    block = parameters.get("bayesian") or {}
    out = {key: parameters[key] for key in BAYESIAN_KEYS if key in parameters}
    out.update({key: value for key, value in block.items() if key != "per_replica"})
    if is_bnn(parameters):
        for key, value in BAYESIAN_DEFAULTS.items():
            out.setdefault(key, value)
        out.setdefault(
            "bayesian_bias", [layer in BAYESIAN_LAYER_TYPES for layer in parameters.get("layer_type")]
        )
        try:
            out.setdefault("kl_warmup_epochs", max(1, min(KL_WARMUP_DEFAULT, int(parameters["epochs"]) // 4)))
        except (KeyError, TypeError, ValueError):  # e.g. a hyperopt search space for epochs
            out.setdefault("kl_warmup_epochs", KL_WARMUP_DEFAULT)
    out["per_replica"] = {**PER_REPLICA_DEFAULTS, **(block.get("per_replica") or {})}
    return out
