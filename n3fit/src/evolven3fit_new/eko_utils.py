from ekobox import gen_theory, gen_op
from eko import run_dglap

from validphys.api import API 
from validphys.loader import Loader
from . import utils


def construct_eko_cards(theoryID, op_card_dict, t_card_dict, q_fin, q_points, x_grid):
    """Return the theory and operator cards used to construct the eko"""
    # theory_card construction
    theory = API.theoryid(theoryid = theoryID).get_description()
    theory.pop("FNS")
    theory.update(t_card_dict)
    t_card = gen_theory.gen_theory_card(theory["PTO"], theory["Q0"], update=theory)
    # construct operator card
    q2_grid = utils.generate_q2grid(
        theory["Q0"],
        q_fin,
        q_points,
        {theory["mb"]: theory["kbThr"], theory["mt"]: theory["ktThr"]},
    )
    op_card = gen_op.gen_op_card(q2_grid, update={"interpolation_xgrid": x_grid})
    op_card.update(op_card_dict)
    return t_card, op_card


def construct_eko_for_fit(t_card, op_card, log, save_path=None):
    """
    Construct the eko operator needed for evolution of fitted pdfs

    Parameters
    ----------
        t_card: dict
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
    eko_op = run_dglap(t_card, op_card)
    if save_path is not None:
        try:
            log.info(f"Saving computed eko to : {save_path}")
            eko_op.dump_tar(save_path)
        except:
            log.error(f"Error saving the eko to : {save_path}")
            pass
    return eko_op
