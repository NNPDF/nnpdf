import logging
from typing import Any, Dict, Optional

import numpy as np

from eko.io import runcards
from eko.matchings import Atlas, nf_default
from eko.quantities.heavy_quarks import MatchingScales

from . import utils

_logger = logging.getLogger(__name__)

EVOLVEN3FIT_CONFIGS_DEFAULTS_TRN = {
    "ev_op_iterations": 1,
    "ev_op_max_order": (1, 0),
    "evolution_method": "truncated",
    "inversion_method": "expanded",
    "n_integration_cores": 1,
}

EVOLVEN3FIT_CONFIGS_DEFAULTS_EXA = {
    "ev_op_max_order": (1, 0),
    "evolution_method": "iterate-exact",
    "inversion_method": "exact",
    "n_integration_cores": 1,
}

NFREF_DEFAULT = 5


def construct_eko_cards(
    nnpdf_theory,
    q_fin,
    q_points,
    x_grid,
    op_card_dict: Optional[Dict[str, Any]] = None,
    theory_card_dict: Optional[Dict[str, Any]] = None,
    legacy40: bool = False,
):
    """
    Return the theory and operator cards used to construct the eko.
    nnpdf_theory is a NNPDF theory card for which we are computing the operator card and eko
    q_fin is the final point of the q grid while q_points is the number of points of the grid.
    x_grid is the x grid to be used.
    op_card_dict and theory_card_dict are optional updates that can be provided respectively to the
    operator card and to the theory card.
    """
    theory, thresholds = load_theory(nnpdf_theory, theory_card_dict)

    mu0 = theory["Q0"]

    # eko needs a value for Qedref and for max nf alphas
    theory["Qedref"] = theory["Qref"]
    theory["MaxNfAs"] = theory["MaxNfPdf"]

    # The Legacy function is able to construct a theory card for eko starting from a NNPDF theory
    theory_card = runcards.Legacy(theory, {}).new_theory

    # construct mugrid

    # Generate the q2grid, if q_fin and q_points are None, use `nf0` to select a default
    q2_grid = utils.generate_q2grid(
        mu0,
        q_fin,
        q_points,
        {
            theory["mc"]: thresholds["c"],
            theory["mb"]: thresholds["b"],
            theory["mt"]: thresholds["t"],
        },
        theory["nf0"],
        legacy40=legacy40,
    )

    masses = np.array([theory["mc"], theory["mb"], theory["mt"]]) ** 2
    thresholds_ratios = np.array([thresholds["c"], thresholds["b"], thresholds["t"]]) ** 2

    atlas = Atlas(
        matching_scales=MatchingScales(masses * thresholds_ratios), origin=(mu0**2, theory["nf0"])
    )

    # Create the eko operator q2grid
    # This is a grid which contains information on (q, nf)
    # in VFNS values at the matching scales need to be doubled so that they are considered in both sides
    ep = 1e-4
    mugrid = []
    for q2 in q2_grid:
        q = float(np.sqrt(q2))
        if nf_default(q2 + ep, atlas) != nf_default(q2 - ep, atlas):
            nf_l = int(nf_default(q2 - ep, atlas))
            nf_u = int(nf_default(q2 + ep, atlas))
            mugrid.append((q, nf_l))
            mugrid.append((q, nf_u))
        else:
            mugrid.append((q, int(nf_default(q2, atlas))))

    # construct operator card
    op_card = build_opcard(op_card_dict, theory, x_grid, mu0, mugrid)

    return theory_card, op_card


def construct_eko_photon_cards(
    nnpdf_theory,
    q_fin,
    x_grid,
    q_gamma,
    op_card_dict: Optional[Dict[str, Any]] = None,
    theory_card_dict: Optional[Dict[str, Any]] = None,
):
    """
    Return the theory and operator cards used to construct the eko_photon.
    nnpdf_theory is a NNPDF theory card for which we are computing the operator card and eko
    q_fin is the final point of the q grid while q_points is the number of points of the grid.
    x_grid is the x grid to be used.
    op_card_dict and theory_card_dict are optional updates that can be provided respectively to the
    operator card and to the theory card.
    """
    theory, _ = load_theory(nnpdf_theory, theory_card_dict)

    # if is eko_photon then mu0 = q_gamma
    mu0 = float(q_gamma)

    # Now make sure the Legacy class still gets a Qedref, which is equal to Qref
    theory["Qedref"] = theory["Qref"]
    theory["MaxNfAs"] = theory["MaxNfPdf"]

    # The Legacy function is able to construct a theory card for eko starting from a NNPDF theory
    theory_card = runcards.Legacy(theory, {}).new_theory

    # The photon needs to be evolved down to Q0
    q_fin = theory["Q0"]
    nf_fin = theory["nf0"]

    # construct mugrid
    mugrid = [(q_fin, nf_fin)]

    # construct operator card
    op_card = build_opcard(op_card_dict, theory, x_grid, mu0, mugrid)

    return theory_card, op_card


def load_theory(nnpdf_theory, theory_card_dict):
    """loads and returns the theory dictionary and the thresholds"""
    if theory_card_dict is None:
        theory_card_dict = {}
    # theory_card construction
    theory = dict(nnpdf_theory)
    theory.pop("FNS")
    theory.update(theory_card_dict)

    if "nfref" not in theory:
        theory["nfref"] = NFREF_DEFAULT

    # Prepare the thresholds according to MaxNfPdf
    thresholds = {"c": theory["kcThr"], "b": theory["kbThr"], "t": theory["ktThr"]}
    if theory["MaxNfPdf"] < 5:
        thresholds["b"] = np.inf
    if theory["MaxNfPdf"] < 6:
        thresholds["t"] = np.inf

    # Setting the thresholds in the theory card to inf if necessary
    theory.update({"kbThr": thresholds["b"], "ktThr": thresholds["t"]})

    return theory, thresholds


def build_opcard(op_card_dict, theory, x_grid, mu0, mugrid):
    """Build the operator card.
    The user provided options should be given as part of ``op_card_dict``
    """
    if op_card_dict is None:
        op_card_dict = {}

    # Taken from cards.py https://github.com/NNPDF/eko/blob/master/src/ekobox/cards.py
    # 7735fdb
    op_card = dict(
        init=(1.65, 4),
        mugrid=[(100.0, 5)],
        xgrid=np.geomspace(1e-7, 1.0, 50).tolist(),
        configs=dict(
            # These three values might be set by op_card_dict
            ev_op_iterations=10,
            n_integration_cores=1,
            polarized=False,
            #
            ev_op_max_order=[10, 0],
            interpolation_polynomial_degree=4,
            interpolation_is_log=True,
            scvar_method=None,
            inversion_method=None,
            evolution_method="iterate-exact",
            time_like=False,
        ),
        debug=dict(skip_singlet=False, skip_non_singlet=False),
    )

    op_card["init"] = (mu0, theory["nf0"])
    op_card["mugrid"] = mugrid
    op_card["xgrid"] = x_grid

    # Specify the evolution options and defaults differently from TRN / EXA
    configs = op_card["configs"]
    if theory.get("ModEv") == "TRN":
        configs.update(EVOLVEN3FIT_CONFIGS_DEFAULTS_TRN)
    elif theory.get("ModEv") == "EXA":
        # Set the default from the theory card unless it was given in the input
        op_card_dict.setdefault("ev_op_iterations", theory.get("IterEv"))

        configs.update(EVOLVEN3FIT_CONFIGS_DEFAULTS_EXA)
        configs["ev_op_iterations"] = op_card_dict["ev_op_iterations"]

    # Note that every entry that is not a dictionary should not be
    # touched by the user and indeed an user cannot touch them
    for key in op_card:
        if key in op_card_dict and isinstance(op_card[key], dict):
            op_card[key].update(op_card_dict[key])
        elif key in op_card_dict:
            _logger.warning("Entry %s is not a dictionary and will be ignored", key)

    op_card = runcards.OperatorCard.from_dict(op_card)
    return op_card
