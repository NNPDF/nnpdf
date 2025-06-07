from collections import defaultdict
import dataclasses
import json
import logging
import pathlib
import sys

from ekobox import apply, genpdf, info_file
import numpy as np

import eko
from eko import basis_rotation, runner
from eko.interpolation import XGrid
from eko.io import manipulate
from validphys.utils import yaml_safe

from . import eko_utils, utils

_logger = logging.getLogger(__name__)

LOG_FILE = "evolven3fit.log"

LOGGING_SETTINGS = {
    "formatter": logging.Formatter("%(asctime)s %(name)s/%(levelname)s: %(message)s"),
    "level": logging.DEBUG,
}


def evolve_fit(
    fit_folder, q_fin, q_points, op_card_dict, theory_card_dict, force, eko_path, dump_eko=None
):
    """
    Evolves all the fitted replica in fit_folder/nnfit

    Parameters
    ----------

        fit_folder: str or pathlib.Path
            path to the folder containing the fit
        q_fin: float
            final point of the q_grid
        q_points: int
            number of points in the q_grid
        op_card_dict: dict
            user settings for the op_card
        theory_card_dict: dict
            user settings for the t_card
        force: bool
            whether to force the evolution to be done again
        eko_path: str or pathlib.Path
            path where the eko is stored (if None the eko will be
            recomputed)
        dump_eko: str or pathlib.Path
            path where the eko is dumped (necessary only if the eko is computed)
    """
    log_file = pathlib.Path(fit_folder) / LOG_FILE
    if log_file.exists():
        if force:
            log_file.unlink()
        else:
            raise FileExistsError(
                f"Log file already exists: {log_file}. Has evolven3fit already been run?"
            )
    log_file = logging.FileHandler(log_file)
    stdout_log = logging.StreamHandler(sys.stdout)
    for log in [log_file, stdout_log]:
        log.setFormatter(LOGGING_SETTINGS["formatter"])

    # The log file will get everything
    log_file.setLevel(LOGGING_SETTINGS["level"])
    # While the terminal only up to info
    stdout_log.setLevel(logging.INFO)

    for logger in (_logger, *[logging.getLogger("eko")]):
        logger.handlers = []
        logger.setLevel(LOGGING_SETTINGS["level"])
        logger.addHandler(log_file)
        logger.addHandler(stdout_log)

    usr_path = pathlib.Path(fit_folder)
    initial_PDFs_dict = load_fit(usr_path)
    x_grid = np.array(initial_PDFs_dict[list(initial_PDFs_dict.keys())[0]]["xgrid"]).astype(float)
    theoryID = utils.get_theoryID_from_runcard(usr_path)

    if eko_path is not None:
        eko_path = pathlib.Path(eko_path)
        _logger.info(f"Loading eko from : {eko_path}")

    if eko_path is None or not eko_path.exists():
        if dump_eko is not None:
            _logger.warning(f"Trying to construct the eko at {dump_eko}")
            theory, op = eko_utils.construct_eko_cards(
                theoryID, q_fin, q_points, x_grid, op_card_dict, theory_card_dict
            )
            runner.solve(theory, op, dump_eko)
            eko_path = dump_eko
        else:
            raise ValueError(f"dump_eko not provided and {eko_path=} not found")

    # Open the EKO in read-only mode, if it needs to be manipulated keep it in memory
    with eko.EKO.read(eko_path) as eko_op:
        # Read the cards directly from the eko to make sure they are consistent
        theory = eko_op.theory_card
        op = eko_op.operator_card
        # And dump them to the log
        _logger.debug(f"Theory card: {json.dumps(theory.raw)}")
        _logger.debug(f"Operator card: {json.dumps(op.raw)}")

        eko_original_xgrid = eko_op.xgrid
        if XGrid(x_grid) != eko_original_xgrid:
            # If the xgrid of the eko is not directly usable, construct a copy in memory
            # by replacing the internal inventory of operators in a readonly copy
            new_xgrid = XGrid(x_grid)
            new_metadata = dataclasses.replace(eko_op.metadata, xgrid=new_xgrid)

            new_operators = {}
            for target_key in eko_op.operators:
                elem = eko_op[target_key.ep]

                if eko_op.metadata.version == "0.13.4":
                    # For eko 0.13.4 xgrid is the internal interpolation so we need to check
                    # whether the rotation is truly needed
                    # <in practice> this means checking whether the operator shape matches the grid
                    oplen = elem.operator.shape[-1]
                    if oplen != len(eko_original_xgrid):
                        # The operator and its xgrid have different shape
                        # either prepare an identity, or this EKO is not supported
                        if oplen != len(x_grid):
                            raise ValueError(
                                f"The operator at {eko_path} is not usable, version not supported"
                            )
                        eko_original_xgrid = XGrid(x_grid)

                new_operators[target_key] = manipulate.xgrid_reshape(
                    elem,
                    eko_original_xgrid,
                    op.configs.interpolation_polynomial_degree,
                    targetgrid=XGrid(x_grid),
                    inputgrid=XGrid(x_grid),
                )

            new_inventory = dataclasses.replace(eko_op.operators, cache=new_operators)
            eko_op = dataclasses.replace(eko_op, metadata=new_metadata, operators=new_inventory)

        # Modify the info file with the fit-specific info
        info = info_file.build(theory, op, 1, info_update={})
        info["NumMembers"] = "REPLACE_NREP"
        info["ErrorType"] = "replicas"
        info["XMin"] = float(x_grid[0])
        info["XMax"] = float(x_grid[-1])
        # Save the PIDs in the info file in the same order as in the evolution
        info["Flavors"] = basis_rotation.flavor_basis_pids
        info.setdefault("NumFlavors", 5)
        dump_info_file(usr_path, info)

        # Read the information from all the sorted replicas into what eko wants
        all_replicas = []
        for rep_dict in initial_PDFs_dict.values():
            # swap photon position to match eko.basis_rotation.flavor_basis_pids
            pdfgrid = np.array(rep_dict["pdfgrid"])
            pdfgrid = np.append(pdfgrid[:, -1].reshape(x_grid.size, 1), pdfgrid[:, :-1], axis=1)
            # and divide by x
            all_replicas.append(pdfgrid.T / x_grid)

        # output is {(Q2, nf): (replica, flavour, x)}
        all_evolved, _ = apply.apply_grids(eko_op, np.array(all_replicas))

        # Now, replica by replica, break into nf blocks
        targetgrid = eko_op.xgrid.tolist()
        by_nf = defaultdict(list)
        for q2, nf in sorted(eko_op.evolgrid, key=lambda ep: ep[1]):
            by_nf[nf].append(q2)
        q2block_per_nf = {nf: sorted(q2s) for nf, q2s in by_nf.items()}

        for replica in initial_PDFs_dict.keys():
            replica_idx = int(replica.split("_")[1]) - 1
            blocks = []
            for nf, q2grid in q2block_per_nf.items():

                def pdf_xq2(pid, x, Q2):
                    x_idx = targetgrid.index(x)
                    pid_idx = info["Flavors"].index(pid)
                    return x * all_evolved[(Q2, nf)][replica_idx][pid_idx][x_idx]

                block = genpdf.generate_block(
                    pdf_xq2,
                    xgrid=targetgrid,
                    sorted_q2grid=q2grid,
                    pids=basis_rotation.flavor_basis_pids,
                )
                blocks.append(block)
            dump_evolved_replica(blocks, usr_path, replica_idx + 1)

    # remove folder:
    # The function dump_evolved_replica uses a temporary folder
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
    for yaml_file in nnfitpath.glob(f"replica_*/{usr_path.name}.exportgrid"):
        data = yaml_safe.load(yaml_file.read_text(encoding="UTF-8"))
        pdf_dict[yaml_file.parent.stem] = data
    return pdf_dict


def dump_evolved_replica(evolved_blocks, usr_path, replica_num):
    """
    Dump the evolved replica given by evolved_block as the replica num "replica_num" in
    the folder usr_path/nnfit/replica_<replica_num>/usr_path.stem.dat

    Parameters
    ----------
        evolved_block: list(numpy.array)
            list of blocks of an evolved PDF
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
        path_where_dump, replica_num, evolved_blocks, pdf_type=to_write_in_head
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
