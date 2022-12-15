from ekobox import gen_theory, gen_op
from eko import run_dglap

from validphys.api import API
from . import utils

from typing import Any, Dict, Optional
import logging

_logger = logging.getLogger(__name__)


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
    theory = API.theoryid(theoryid=theoryID).get_description()
    theory.pop("FNS")
    theory.update(theory_card_dict)
    theory_card = gen_theory.gen_theory_card(theory["PTO"], theory["Q0"], update=theory)
    # construct operator card
    q2_grid = utils.generate_q2grid(
        theory["Q0"],
        q_fin,
        q_points,
        {theory["mb"]: theory["kbThr"], theory["mt"]: theory["ktThr"]},
    )
    op_card = gen_op.gen_op_card(q2_grid, update={"interpolation_xgrid": x_grid})
    op_card.update(op_card_dict)
    return theory_card, op_card


def construct_eko_for_fit(theory_card, op_card, save_path=None):
    """
    Construct the eko operator needed for evolution of fitted pdfs

    Parameters
    ----------
        theory_card: dict
            theory card to use for the eko
        op_card: dict
            operator card to use for the eko
        save_path: pathlib.Path
            path where the eko will be saved (the eko
            won't be saved if save_path is None)
    Returns
    -------
        : eko.output.Output
        eko operator
    """
    # generate eko operator
    if save_path is not None:
        if not save_path.parent.exists():
            raise FileNotFoundError(
                f"Path where eko should be dumped does not exist: {save_path}"
            )
    eko_op = run_dglap(theory_card, op_card)
    if save_path is not None:
        # Here we want to catch all possible exceptions in order to avoid losing the computed eko
        try:
            _logger.info(f"Saving computed eko to : {save_path}")
            eko_op.dump_tar(save_path)
        except:
            _logger.error(f"Error saving the eko to : {save_path}")
            pass
    return eko_op
