import pathlib
import logging
import sys

import numpy as np
from reportengine.compat import yaml

from ekobox import genpdf, gen_info
from ekomark import apply
from eko import basis_rotation as br
from eko import output
from validphys.api import API

from . import utils, eko_utils


log = logging.getLogger(__name__)

LOG_FILE = "evolven3fit_new.log"

LOGGING_SETTINGS = {"formatter" : "%(asctime)s %(name)s/%(levelname)s: %(message)s",
                    "level" : logging.INFO
}

def evolve_fit(
    conf_folder,
    q_fin,
    q_points,
    op_card_dict,
    t_card_dict,
    force,
    eko_path=None,
    dump_eko=None,
):
    """
    Evolves all the fitted replica in conf_folder/nnfit

    Parameters
    ----------

        conf_folder: str or pathlib.Path
            path to the folder containing the fit
        q_fin: float
            final point of the q_grid
        q_points: int
            number of points in the q_grid
        op_card_dict: dict
            user settings for the op_card
        t_card_dict: dict
            user settings for the t_card
        force: bool
            whether to force the evolution to be done again
        eko_path: str or pathlib.Path
            path where the eko is stored (if None the eko will be
            recomputed)
        dump_eko: str or pathlib.Path
            path where the eko is dumped (if None the eko won't be
            stored)
    """
    log_file = pathlib.Path(conf_folder) / LOG_FILE
    if log_file.exists():
        if force:
            log_file.unlink()
        else:
            raise FileExistsError(
                f"Log file already exists: {log_file} evolven3fit_new has already been run?"
            )
    log_file = logging.FileHandler(log_file)
    stdout_log = logging.StreamHandler(sys.stdout)
    for log in [log_file, stdout_log]:
        log.setLevel(LOGGING_SETTINGS["level"])
        log.setFormatter(
            logging.Formatter(LOGGING_SETTINGS["formatter"])
        )
    for logger in (log, *[logging.getLogger("eko")]):
        logger.handlers = []
        logger.setLevel(LOGGING_SETTINGS["level"])
        logger.addHandler(log_file)
        logger.addHandler(stdout_log)

    usr_path = pathlib.Path(conf_folder)
    initial_PDFs_dict = load_fit(usr_path)
    x_grid = np.array(
        initial_PDFs_dict[list(initial_PDFs_dict.keys())[0]]["xgrid"]
    ).astype(np.float)
    theoryID = utils.get_theoryID_from_runcard(usr_path)
    theory, op = eko_utils.construct_eko_cards(
        theoryID, q_fin, q_points, x_grid, op_card_dict, t_card_dict
    )
    if eko_path is not None:
        eko_path = pathlib.Path(eko_path)
        log.info(f"Loading eko from : {eko_path}")
        eko_op = output.Output.load_tar(eko_path)
    else:
        try:
            log.info(f"Loading eko from theory {theoryID}")
            theory_eko_path = (API.theoryid(theoryid = theoryID).path)/'eko.tar'
            eko_op = output.Output.load_tar(theory_eko_path)
        except:
            log.info(f"eko not found in theory {theoryID}, we will construct it")
            eko_op = eko_utils.construct_eko_for_fit(theory, op, log, dump_eko)
            pass   
    eko_op.xgrid_reshape(targetgrid=x_grid, inputgrid=x_grid)
    info = gen_info.create_info_file(theory, op, 1, info_update={})
    info["NumMembers"] = "REPLACE_NREP"
    info["ErrorType"] = "replicas"
    info["AlphaS_Qs"] = list(eko_op["Q2grid"])
    info["XMin"] = float(x_grid[0])
    info["XMax"] = float(x_grid[-1])
    dump_info_file(usr_path, info)
    for replica in initial_PDFs_dict.keys():
        evolved_block = evolve_exportgrid(initial_PDFs_dict[replica], eko_op, x_grid)
        dump_evolved_replica(
            evolved_block, usr_path, int(replica.removeprefix("replica_"))
        )
    # remove folder:
    # The function dump_evolved_replica dumps the replica files in a temporary folder
    # We need then to remove it after fixing the position of those replica files
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
    pdf_dict = {}
    for yaml_file in nnfitpath.glob("replica_*/*.exportgrid"):
        data = yaml.safe_load(yaml_file.read_text(encoding="UTF-8"))
        pdf_dict[yaml_file.parent.stem] = data
    return pdf_dict


def evolve_exportgrid(exportgrid, eko, x_grid):
    """
    Evolves the provided exportgrid for the desired replica with the eko and returns the evolved block

    Parameters
    ----------
        exportgrid: dict
            exportgrid of pdf at fitting scale
        eko: eko object
            eko operator for evolution
        xgrid: list
            xgrid to be used as the targetgrid
    Returns
    -------
        : np.array
        evolved block
    """
    # construct LhapdfLike object
    pdf_grid = np.array(exportgrid["pdfgrid"]).transpose()
    pdf_to_evolve = utils.LhapdfLike(pdf_grid, exportgrid["q20"], x_grid)
    # evolve pdf
    evolved_pdf = apply.apply_pdf(eko, pdf_to_evolve)
    # generate block to dump
    targetgrid = list(eko["targetgrid"])

    def ev_pdf(pid, x, Q2):
        return x * evolved_pdf[Q2]["pdfs"][pid][targetgrid.index(x)]

    block = genpdf.generate_block(
        ev_pdf,
        xgrid=targetgrid,
        Q2grid=list(eko["Q2grid"]),
        pids=br.flavor_basis_pids,
    )
    return block


def dump_evolved_replica(evolved_block, usr_path, replica_num):
    """
    Dump the evolved replica given by evolved_block as the replica num "replica_num" in
    the folder usr_path/nnfit/replica_<replica_num>/usr_path.stem.dat

    Parameters
    ----------
        evolved_block: numpy.array
            block of an evolved PDF
        usr_path: pathlib.Path
            path of the fit folder
        replica_num: int
            replica number
    """
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    # create folder to dump the evolved replica if it does not exist
    path_where_dump.mkdir(exist_ok=True)
    to_write_in_head = f"PdfType: replica\nFromMCReplica: {replica_num}\n"
    genpdf.export.dump_blocks(
        path_where_dump, replica_num, [evolved_block], pdf_type=to_write_in_head
    )
    # fixing_replica_path
    utils.fix_replica_path(usr_path, replica_num)


def dump_info_file(usr_path, info):
    """
    Dump the info file given by info in the folder usr_path/nnfit/usr_path.stem.info.

    Parameters
    ----------
        usr_path: pathlib.Path
            path of the fit folder
        info: dict
            info of the fit
    """
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    genpdf.export.dump_info(path_where_dump, info)
    utils.fix_info_path(usr_path)
