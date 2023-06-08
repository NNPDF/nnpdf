import logging
from typing import Any, Dict, Optional

from ekobox.cards import _operator as default_op_card
import numpy as np

from eko.io import runcards
from eko.matchings import Atlas, nf_default
from eko.quantities.heavy_quarks import MatchingScales
from validphys.loader import Loader

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
    "ev_op_iterations": 10,
    "ev_op_max_order": (1, 0),
    "evolution_method": "iterate-exact",
    "inversion_method": "exact",
    "n_integration_cores": 1,
}

NFREF_DEFAULT = 5
NF0_DEFAULT = 4


def construct_eko_cards(
    theoryID,
    q_fin,
    q_points,
    x_grid,
    op_card_dict: Optional[Dict[str, Any]] = None,
    theory_card_dict: Optional[Dict[str, Any]] = None,
):
    """
    Return the theory and operator cards used to construct the eko.
    theoryID is the ID of the theory for which we are computing the theory and operator card.
    q_fin is the final point of the q grid while q_points is the number of points of the grid.
    x_grid is the x grid to be used.
    op_card_dict and theory_card_dict are optional updates that can be provided respectively to the
    operator card and to the theory card.
    """
    if theory_card_dict is None:
        theory_card_dict = {}
    if op_card_dict is None:
        op_card_dict = {}

    # theory_card construction
    theory = Loader().check_theoryID(theoryID).get_description()
    theory.pop("FNS")
    theory.update(theory_card_dict)
    if "nfref" not in theory:
        theory["nfref"] = NFREF_DEFAULT
    if "nf0" not in theory:
        theory["nf0"] = NF0_DEFAULT
        
    # Prepare the thresholds according to MaxNfPdf
    thresholds = {"c": theory["kcThr"], "b": theory["kbThr"], "t": theory["ktThr"]}
    if theory["MaxNfPdf"] < 5:
        thresholds["b"] = np.inf
    if theory["MaxNfPdf"] < 6:
        thresholds["t"] = np.inf
        
    # Setting the thresholds in the theory card to inf if necessary
    theory.update({"kbThr":thresholds["b"], "ktThr":thresholds["t"] })
    
    # The Legacy function is able to construct a theory card for eko starting from an NNPDF theory
    legacy_class = runcards.Legacy(theory, {})
    theory_card = legacy_class.new_theory

    # construct operator card
    q2_grid = utils.generate_q2grid(
        theory["Q0"],
        q_fin,
        q_points,
        {theory["mb"]: thresholds["b"], theory["mt"]: thresholds["t"]},
    )
    op_card = default_op_card
    masses = np.array([theory["mc"], theory["mb"], theory["mt"]]) ** 2
    thresholds_ratios = np.array([thresholds["c"], thresholds["b"], thresholds["t"]]) ** 2

    atlas = Atlas(
        matching_scales=MatchingScales(masses * thresholds_ratios),
        origin=(theory["Q0"] ** 2, theory["nf0"]),
    )
    op_card.update(
        {
            "mu0": theory["Q0"],
            "mugrid": [(float(np.sqrt(q2)), int(nf_default(q2, atlas))) for q2 in q2_grid],
        }
    )
    op_card["xgrid"] = x_grid
    # Specific defaults for evolven3fit evolution
    if theory["ModEv"] == "TRN":
        op_card["configs"].update(EVOLVEN3FIT_CONFIGS_DEFAULTS_TRN)
    if theory["ModEv"] == "EXA":
        op_card["configs"].update(EVOLVEN3FIT_CONFIGS_DEFAULTS_EXA)
    # User can still change the configs via op_card_dict

    # Note that every entry that is not a dictionary should not be
    # touched by the user and indeed an user cannot touch them
    for key in op_card:
        if key in op_card_dict and isinstance(op_card[key], dict):
            op_card[key].update(op_card_dict[key])
        elif key in op_card_dict:
            _logger.warning("Entry %s is not a dictionary and will be ignored", key)

    op_card = runcards.OperatorCard.from_dict(op_card)
    return theory_card, op_card


def split_evolgrid(evolgrid):
    """Split the evolgrid in blocks according to the number of flavors and repeating the last entry of one block in the first entry of the next."""
    evolgrid_index_list = []
    evolgrid.sort()
    starting_nf = evolgrid[0][1]
    for evo_point in evolgrid:
        current_nf = evo_point[1]
        if current_nf != starting_nf:
            evolgrid_index_list.append(evolgrid.index(evo_point))
            starting_nf = current_nf
    start_index = 0
    evolgrid_list = []
    for index in evolgrid_index_list:
        evolgrid_list.append(evolgrid[start_index : index + 1])
        start_index = index
    evolgrid_list.append(evolgrid[start_index:])
    return evolgrid_list
