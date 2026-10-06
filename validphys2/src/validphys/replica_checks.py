"""
replica_checks.py

Per-replica checks for fits where several replicas come out of a single training,
i.e. the pseudo-replicas sampled from a Bayesian neural network. In that case the
positivity status that n3fit's stopping object computes during training belongs to
the trained model only, not to the individual posterior samples written to
``nnfit/replica_*``. The functions here evaluate positivity for each replica
separately, with the same definition n3fit uses during training:

    loss_set = lambda_set * sum_i elu(-prediction_i, alpha=1e-7)
    replica passes  <=>  loss_set < threshold for every positivity set

(see :py:class:`n3fit.layers.losses.LossPositivity` and
:py:class:`n3fit.stopping.Positivity`). n3fit's stopping evaluates this on the
*validation* model, whose positivity layers are separate instances from the training
model's (``ObservableWrapper.__call__`` builds a new ``LossPositivity`` per model) and are
never updated by ``LagrangeCallback``: their multiplier stays at
``parameters::positivity::initial`` for the whole fit, while only the training penalty
ramps towards ``maxlambda``. :py:func:`stopping_multipliers` reproduces that, so a sample
is held to exactly the check NNPDF applies to a replica.

Integrability is not handled here: the integrability numbers written to each
replica's ``.json`` are already computed from that replica's own PDF and are
vetoed by postfit (see :py:mod:`validphys.fitveto`).

Runcard switches (``parameters::bayesian::per_replica``, see ``n3fit.bayesian_settings``):
``positivity: true`` makes n3fit write a per-sample ``pos_state`` (postfit's standard
convergence veto then acts per sample); ``postfit_veto: true`` makes postfit add an explicit
per-replica positivity veto, using n3fit's per-sample result where available and
recomputing positivity here from the ``.exportgrid`` otherwise. ``postfit --per-replica-checks``
does the latter for fits run without the switch.
"""

import logging
import pathlib

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from validphys.convolution import fk_predictions
from validphys.core import PDF
from validphys.utils import yaml_safe

log = logging.getLogger(__name__)

# n3fit defaults: alpha of LossPositivity, the positivity threshold of the stopping object
# and how often the positivity multipliers are updated (n3fit.model_trainer.PUSH_POSITIVITY_EACH;
# duplicated here so that validphys does not import n3fit)
POSITIVITY_ALPHA = 1e-7
POSITIVITY_THRESHOLD = 1e-6
PUSH_POSITIVITY_EACH = 100

# key written to each replica's .json when its pos_state was computed per replica
PER_REPLICA_JSON_KEY = "per_replica_positivity"
# runcard switch for the postfit veto: parameters::bayesian::per_replica::postfit_veto


def postfit_veto_requested(runcard):
    """True if the fit runcard asks postfit for the per-replica positivity veto"""
    params = runcard.get("parameters") or {}
    per_replica = (params.get("bayesian") or {}).get("per_replica") or {}
    return bool(per_replica.get("postfit_veto", False))

EXPORTGRID_PDG = {
    "TBAR": -6,
    "BBAR": -5,
    "CBAR": -4,
    "SBAR": -3,
    "UBAR": -2,
    "DBAR": -1,
    "GLUON": 21,
    "D": 1,
    "U": 2,
    "S": 3,
    "C": 4,
    "B": 5,
    "T": 6,
    "PHT": 22,
}


def _elu(x, alpha):
    return np.where(x > 0, x, alpha * np.expm1(np.minimum(x, 0)))


def _lm_initial_and_multiplier(initial, multiplier, max_lambda, steps):
    """Same as :py:func:`n3fit.model_trainer._LM_initial_and_multiplier`"""
    if multiplier is None:
        if initial is None:
            initial = 1.0
        multiplier = pow(max_lambda / initial, 1 / max(steps, 1))
    elif initial is None:
        initial = max_lambda / pow(multiplier, steps)
    return initial, multiplier


def stopping_multipliers(pos_dicts, parameters):
    """Positivity multiplier of each set as used by n3fit's stopping positivity check.

    That check reads the validation model's losses, whose ``LossPositivity`` layers keep
    the initial multiplier (``parameters::positivity::initial``, or the value derived from
    ``multiplier``/``maxlambda`` as in ``n3fit.model_trainer._LM_initial_and_multiplier``)
    for the whole fit: ``LagrangeCallback`` only ramps the training model's layers.
    """
    pos_par = parameters.get("positivity") or {}
    steps = int(parameters["epochs"] / PUSH_POSITIVITY_EACH)
    lambdas = {}
    for pos in pos_dicts:
        initial, _ = _lm_initial_and_multiplier(
            pos_par.get("initial"), pos_par.get("multiplier"), pos["lambda"], steps
        )
        lambdas[pos["name"]] = initial
    return pd.Series(lambdas)


def positivity_losses(pdf, pos_dicts, lambdas=None, alpha=POSITIVITY_ALPHA):
    """Positivity loss of every PDF member for every positivity set.

    Parameters
    ----------
    pdf: validphys.core.PDF
        any PDF accepted by :py:func:`validphys.convolution.fk_predictions` whose
        members are the replicas to check (e.g. ``n3fit.vpinterface.N3PDF`` or
        :py:class:`ExportGridPDF`)
    pos_dicts: list(dict)
        positivity sets as produced by
        :py:func:`validphys.n3fit_data.posdatasets_fitting_pos_dict`
    lambdas: pd.Series or pd.DataFrame
        multiplier per set (Series), or per member and set (DataFrame with one row per
        member), see :py:func:`stopping_multipliers`. Defaults to each set's ``maxlambda``.

    Returns
    -------
    pd.DataFrame
        losses with one row per PDF member and one column per positivity set
    """
    losses = {}
    for pos in pos_dicts:
        fktables = [fk for ds in pos["datasets"] for fk in ds.fktables_data]
        if len(fktables) != 1:
            raise NotImplementedError(
                f"Positivity set {pos['name']} has {len(fktables)} FK tables, expected 1"
            )
        preds = fk_predictions(fktables[0], pdf).to_numpy()  # (ndata, nmembers)
        losses[pos["name"]] = _elu(-preds, alpha).sum(axis=0)
    losses = pd.DataFrame(losses)
    if lambdas is None:
        lambdas = pd.Series({pos["name"]: pos["lambda"] for pos in pos_dicts})
    if isinstance(lambdas, pd.DataFrame):
        return losses * lambdas[losses.columns].to_numpy()
    return losses * lambdas[losses.columns]


def positivity_states(losses, threshold=POSITIVITY_THRESHOLD):
    """Boolean array, True for the members passing every positivity set"""
    return (losses < threshold).all(axis=1).to_numpy()


class ExportGridSet:
    """LHAPDF-like set built from n3fit ``.exportgrid`` files (PDF at the fitting
    scale), interpolated in log(x) with a monotone (PCHIP) interpolant, which
    cannot overshoot below zero where the PDFs vanish. Only the ``grid_values``
    interface used by :py:func:`validphys.convolution.fk_predictions` is provided."""

    def __init__(self, exportgrids):
        self.q0 = exportgrids[0]["q20"] ** 0.5
        self.xgrid = np.asarray(exportgrids[0]["xgrid"])
        labels = exportgrids[0]["labels"]
        self.pdg = [EXPORTGRID_PDG[l] for l in labels]
        # (members, x, flavours)
        self.values = np.array([eg["pdfgrid"] for eg in exportgrids])
        self.splines = PchipInterpolator(np.log(self.xgrid), self.values, axis=1)

    def grid_values(self, flavors, xgrid, qgrid):
        qgrid = np.atleast_1d(qgrid)
        if not np.allclose(qgrid, self.q0, rtol=1e-3):
            log.warning(
                "Exportgrids are at Q0=%s, but were queried at Q=%s; using Q0", self.q0, qgrid
            )
        xgrid = np.atleast_1d(xgrid)
        vals = self.splines(np.log(np.clip(xgrid, self.xgrid[0], self.xgrid[-1])))
        idx = [self.pdg.index(int(f)) if int(f) in self.pdg else None for f in flavors]
        out = np.zeros((self.values.shape[0], len(flavors), len(xgrid)))
        for i, j in enumerate(idx):
            if j is not None:
                out[:, i, :] = vals[:, :, j]
        return np.repeat(out[..., np.newaxis], len(qgrid), axis=-1)


class ExportGridPDF(PDF):
    """validphys PDF whose members are the replicas behind a list of exportgrids"""

    def __init__(self, exportgrids, name="exportgrid_replicas"):
        super().__init__(name)
        self._set = ExportGridSet(exportgrids)
        self._info = {"ErrorType": "replicas", "NumMembers": len(exportgrids)}

    def load(self):
        return self._set

    def load_t0(self):
        return self._set

    def get_members(self):
        return self._info["NumMembers"]


def load_exportgrid(replica_path, fitname):
    path = pathlib.Path(replica_path) / f"{fitname}.exportgrid"
    return yaml_safe.load(path.read_text(encoding="utf-8"))


def runcard_pos_dicts(runcard):
    """Positivity sets (with their FK tables) and threshold of a fit runcard, in the
    format of :py:func:`validphys.n3fit_data.posdatasets_fitting_pos_dict`"""
    from validphys.api import API
    from validphys.loader import Loader
    from validphys.n3fit_data import _fitting_lagrange_dict

    loader = Loader()
    theoryid = runcard["theory"]["theoryid"]
    # same cuts as in the fit: internal rules plus the runcard's filter-rule settings
    # (e.g. NNPDF4.0 imposes PDF positivity only for x > 0.1 via added_filter_rules)
    rule_keys = (
        "datacuts",
        "added_filter_rules",
        "filter_rules",
        "drop_internal_rules",
        "default_filter_settings",
        "default_filter_rules",
    )
    rules = API.rules(
        theoryid=theoryid, use_cuts="internal", **{k: runcard[k] for k in rule_keys if k in runcard}
    )
    pos_dicts = []
    for posset in runcard["positivity"]["posdatasets"]:
        spec = loader.check_posset(theoryid, posset["dataset"], float(posset["maxlambda"]), rules)
        pos_dicts.append(_fitting_lagrange_dict(spec))
    if not pos_dicts:
        raise ValueError("No positivity datasets found in the fit runcard")
    threshold = runcard.get("parameters", {}).get("positivity", {}).get("threshold")
    return pos_dicts, (POSITIVITY_THRESHOLD if threshold is None else threshold)


def replicas_positivity(replica_paths, fitname, runcard, chunk_size=200):
    """Recompute per-replica positivity from the ``.exportgrid`` of each replica, with
    the multipliers of n3fit's stopping check (:py:func:`stopping_multipliers`).

    Returns a boolean array (True = passes) and the losses DataFrame, one row per path.
    """
    pos_dicts, threshold = runcard_pos_dicts(runcard)
    lambdas = stopping_multipliers(pos_dicts, runcard["parameters"])
    all_losses = []
    for start in range(0, len(replica_paths), chunk_size):
        chunk = replica_paths[start : start + chunk_size]
        pdf = ExportGridPDF([load_exportgrid(p, fitname) for p in chunk])
        all_losses.append(positivity_losses(pdf, pos_dicts, lambdas))
    losses = pd.concat(all_losses, ignore_index=True)
    return positivity_states(losses, threshold), losses
