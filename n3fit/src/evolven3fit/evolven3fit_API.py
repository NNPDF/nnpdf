"""
The evolven3fit_API module provides the necessary tools to evolve a fit 
using eko and to prepare it for postfit.  
"""

import pathlib
import logging
import numpy as np
import yaml
import os
import glob

from ekobox import gen_theory, gen_op, genpdf, gen_info
from ekomark import apply
from eko import basis_rotation as br
from eko import run_dglap
from eko import output

from . import utils
from validphys.loader import Loader


logger = logging.getLogger(__name__)


def evolve_fit(conf_folder, op_card_dict, t_card_dict, eko_path=None, dump_eko=None):
    """
    Evolves all the fitted replica in conf_folder/nnfit

    Parameters
    ----------

        conf_folder: str or pathlib.Path
            path to the folder containing the fit
        op_card_dict: dict
            user settings for the op_card
        t_card_dict: dict
            user settings for the t_card
        eko_path: str or pathlib.Path
            path where the eko is stored (if None the eko will be
            recomputed)
        dump_eko: str or pathlib.Path
            path where the eko is dumped (if None the eko won't be
            stored)
    """
    log_file = pathlib.Path(conf_folder) / "evolven3fit.log"
    if log_file.exists():
        raise SystemError("Log file already exists: evolven3fit has already been run?")
    log_file = logging.FileHandler(log_file)
    log_file.setLevel(logging.INFO)
    log_file.setFormatter(
        logging.Formatter("%(asctime)s %(name)s/%(levelname)s: %(message)s")
    )
    for logger_ in (logger, *[logging.getLogger("eko")]):
        logger_.handlers = []
        logger_.setLevel(logging.INFO)
        logger_.addHandler(log_file)
    usr_path = pathlib.Path(conf_folder)
    initial_PDFs_dict = load_fit(usr_path)
    theory, op = construct_eko_cards(usr_path, op_card_dict, t_card_dict)
    if eko_path is not None:
        eko_path = pathlib.Path(eko_path)
        eko_op = output.Output.load_tar(eko_path)
    else:
        eko_op = construct_eko_for_fit(theory, op, dump_eko)
    info = gen_info.create_info_file(theory, op, 1, info_update={})
    info["NumMembers"] = "REPLACE_NREP"
    info["ErrorType"] = "replicas"
    dump_info_file(usr_path, info)
    utils.fix_info_path(usr_path)
    for replica in initial_PDFs_dict.keys():
        evolved_block = evolve_exportgrid(
            initial_PDFs_dict[replica], eko_op, theory, op
        )
        dump_evolved_replica(
            evolved_block, usr_path, int(replica.removeprefix("replica_"))
        )
        # fixing_replica_path
        utils.fix_replica_path(usr_path, int(replica.removeprefix("replica_")))
    # remove folder
    (usr_path / "nnfit" / usr_path.stem).rmdir()


def load_fit(usr_path):
    """
    Loads all the replica pdfs at fitting scale in usr_path and return the exportgrids

    Parameters
    ----------

        usr_path: pathlib.Path
            path to the folder containing the fit

    Returns
    -------

            : dict
            exportgrids info
    """
    nnfitpath = usr_path / "nnfit"
    replica_list = []
    for file in nnfitpath.glob("replica_*"):
            replica_list.append(file)
    pdf_dict = {}
    for replica in replica_list:
        yaml_file = pathlib.Path(replica) / (usr_path.stem + ".exportgrid")
        with yaml_file.open() as fp:
            data = yaml.safe_load(fp)
        pdf_dict[replica.stem] = data
    return pdf_dict


def construct_eko_cards(usr_path, op_card_dict, t_card_dict):
    """Return the theory and operator cards used to construct the eko"""
    # read the runcard
    my_runcard = utils.read_runcard(usr_path)
    # theory_card construction
    theory = Loader().check_theoryID(my_runcard["theory"]["theoryid"]).get_description()
    theory.pop("FNS")
    theory.update(t_card_dict)
    t_card = gen_theory.gen_theory_card(theory["PTO"], theory["Q0"], update=theory)
    # construct operator card
    op_x_grid = utils.generate_x_grid()
    q2_grid = utils.generate_q2grid(theory["Q0"], 1.0e5)
    op_card = gen_op.gen_op_card(q2_grid, update={"interpolation_xgrid": op_x_grid})
    op_card.update(op_card_dict)
    return t_card, op_card


# Temporary solution. Then it will be loaded from the theory itself
def construct_eko_for_fit(t_card, op_card, save_path=None):
    """
    Construct the eko operator needed for evolution of fitted pdfs

    Parameters
    ----------
        usr_path: pathlib.Path
            path to the folder containing the fit
        save_path: pathlib.Path
            path where the eko will be saved (the eko
            won't be saved if save_path is None)
    Returns
    -------
        : eko.output.Output
        eko operator
        : dict
        theory runcard
        : dict
        operator runcard
    """

    # generate eko operator (temporary because it will be loaded from theory)
    eko_op = run_dglap(t_card, op_card)
    if save_path is not None:
        try:
            eko_op.dump_tar(save_path)
        except:
            pass
    return eko_op


def evolve_exportgrid(exportgrid, eko, theory_card, operator_card):
    """
    Evolves the provided exportgrid for the desired replica with the eko and returns the evolved block

    Parameters
    ----------
        exportgrid: dict
            exportgrid of pdf at fitting scale
        eko: eko object
            eko operator for evolution
        theory_card: dict
            theory card
        operator_card: dict
            operator card
    Returns
    -------
        : np.array
        evolved block
    """
    # construct LhapdfLike object
    pdf_grid = np.array(exportgrid["pdfgrid"]).transpose()
    x_grid = np.array(exportgrid["xgrid"]).astype(np.float)
    pdf_to_evolve = utils.LhapdfLike(pdf_grid, exportgrid["q20"], x_grid)
    # evolve pdf
    evolved_pdf = apply.apply_pdf(eko, pdf_to_evolve)
    # generate block to dump
    targetgrid = operator_card["interpolation_xgrid"]
    block = genpdf.generate_block(
        lambda pid, x, Q2, evolved_PDF=evolved_pdf: x
        * evolved_PDF[Q2]["pdfs"][pid][targetgrid.index(x)],
        xgrid=targetgrid,
        Q2grid=operator_card["Q2grid"],
        pids=br.flavor_basis_pids,
    )
    return block


def dump_evolved_replica(evolved_block, usr_path, replica_num):
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    # create folder to dump the evolved replica if it does not exist
    if not os.path.exists(path_where_dump):
        os.makedirs(path_where_dump)
    to_write_in_head = "PdfType: replica\nFromMCReplica: " + str(replica_num) + "\n"
    genpdf.export.dump_blocks(
        path_where_dump, replica_num, [evolved_block], pdf_type=to_write_in_head
    )


def dump_info_file(usr_path, info):
    # dump info file if not already there
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    info_path = path_where_dump / (usr_path.stem + ".info")
    if not info_path.is_file():
        genpdf.export.dump_info(path_where_dump, info)
