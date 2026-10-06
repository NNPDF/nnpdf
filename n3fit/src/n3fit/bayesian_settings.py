"""
bayesian_settings.py

All options of a Bayesian-NN fit live in one block of the runcard, ``parameters::bayesian``
(the architecture itself, ``parameters::layer_type`` with ``VBDense`` layers, stays where it is):

    parameters:
      layer_type: [dense, dense, VBDense]
      bayesian:
        n_bnn_models: 1             # trained BNN models (one per SLURM task)
        n_bnn_samples: 1200         # posterior samples written per model
        sampling_space: weight      # weight | function (linearized Laplace)
        prior_prec: 1.0
        bayesian_bias: [false, false, true]
        bayesian_preproc: false
        bayes_dropout: 0.0
        bayesian_flow: false        # plus flow_n_couplings, flow_hidden_units
        rank: 4                     # VBDense_correlated only, plus u_init
        per_replica:                # per-sample quantities instead of copies of the trained model's
          chi2: false               # chi2, erf_tr, erf_vl of each sample (postfit's chi2 veto then acts per sample)
          positivity: false         # pos_state of each sample, with n3fit's own stopping positivity check
          postfit_veto: false       # postfit: explicit per-replica positivity veto (recomputed from the
                                    #   exportgrid for replicas that were not checked by n3fit)
          posterior_stopping_samples: 0   # > 0: stopping requires positivity at the posterior mean
                                          #   and at this many fixed posterior draws
          rejection_sampling: false       # keep only posterior draws that pass positivity
          rejection_max_draws: null       # default 50 * n_bnn_samples

Runcards that define these keys directly under ``parameters`` (the format used so far) keep
working: those values are read when the ``bayesian`` block does not define them.

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
)
PER_REPLICA_DEFAULTS = {
    "chi2": False,
    "positivity": False,
    "postfit_veto": False,
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


def bayesian_parameters(parameters):
    """
    The Bayesian options of a fit, from ``parameters::bayesian`` with a fallback to the same
    keys directly under ``parameters``. Keys that are set nowhere are absent, so that callers
    keep their own defaults (``bayesian_parameters(p).get("prior_prec", ...)``).
    ``per_replica`` is always present, filled with ``PER_REPLICA_DEFAULTS``.
    """
    block = parameters.get("bayesian") or {}
    out = {key: parameters[key] for key in BAYESIAN_KEYS if key in parameters}
    out.update({key: value for key, value in block.items() if key != "per_replica"})
    out["per_replica"] = {**PER_REPLICA_DEFAULTS, **(block.get("per_replica") or {})}
    return out
