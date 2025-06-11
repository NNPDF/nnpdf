from collections import defaultdict
import dataclasses
import json
import logging
import pathlib
import shutil
import sys
import tempfile

from ekobox import apply, genpdf, info_file
import numpy as np

import eko
from eko import basis_rotation, runner
from eko.interpolation import XGrid
from eko.io import manipulate
from validphys.pdfbases import PIDS_DICT
from validphys.utils import yaml_safe

_logger = logging.getLogger(__name__)

LOG_FILE = "evolven3fit.log"

LOGGING_SETTINGS = {
    "formatter": logging.Formatter("%(asctime)s %(name)s/%(levelname)s: %(message)s"),
    "level": logging.DEBUG,
}

LABEL_TO_PIDS = {v: k for k, v in PIDS_DICT.items()}


@dataclasses.dataclass
class ExportGrid:
    """
    Holds information about the PDF at a fixed scale.

    q20: float
        Value of Q2 at which the PDF is evaluated

    xgrid: np.ndarray
        The points in x at which the PDF is evaluated

    pdfgrid: np.ndarray
        A list of the x*PDF values for all flavours for every point in x
        shape (len(xgrid), len(labels))

    labels: list[str]
        A list of the flavours contained in each element of the pdfgrid,
        defaults to: ['TBAR', 'BBAR', 'CBAR', 'SBAR', 'UBAR', 'DBAR', 'GLUON', 'D', 'U', 'S', 'C', 'B', 'T', 'PHT']

    replica: int
        Index of the corresponding monte carlo replica, should be ``None`` for hessian fits
    """

    q20: float
    xgrid: np.ndarray
    pdfgrid: np.ndarray
    labels: list = (
        'TBAR',
        'BBAR',
        'CBAR',
        'SBAR',
        'UBAR',
        'DBAR',
        'GLUON',
        'D',
        'U',
        'S',
        'C',
        'B',
        'T',
        'PHT',
    )
    replica: int = None

    def __post_init__(self):
        """Convert possible lists to arrays"""
        self.pdfgrid = np.array(self.pdfgrid)
        self.xgrid = np.array(self.xgrid)

    @property
    def pids(self):
        """Return PIDs instead of labels"""
        return [LABEL_TO_PIDS[i] for i in self.labels]

    @property
    def hessian(self):
        """If self.replica is None, assume the fit is hessian"""
        return self.replica is None

    @property
    def pdfvalues(self):
        """Return the PDF, i.e., pdfgrid / xgrid,
        with a shape (flavours, xgrid)
        """
        return self.pdfgrid.T / self.xgrid


def evolve_exportgrid(eko_path, exportgrids):
    """Takes the path to an EKO and a list of exportgrids,
    returns a tuple with an info file and the
    evolved exportgrid as a dictionary of the form:

        {
            (Q_1^2, nf1): (replica, flavours, x),
            (Q_2^2, nf1): (replica, flavours, x),
            ...
            (Q_3^2, nf2): (replica, flavours, x),
        }

        with the output grouped by nf and sorted in ascending order by Q2

    Parameters:
        eko_path: pathlib.Path
            Path to the evolution eko
        exportgrids: list[ExportGrid]

    Returns
        info_file: eko_box.info_file
            dict-like object with the info file information

        evolved_replicas: dict

    """
    # Check that all exportgrid objects have been evaluated for 1) The same value of Q, the same value of x
    ref = exportgrids[0]

    hessian_fit = ref.replica is None
    for egrid in exportgrids:
        assert egrid.q20 == ref.q20, "Different values of q0 found among the exportgrids"
        np.testing.assert_allclose(
            ref.xgrid, egrid.xgrid, err_msg="ExportGrids are not all evaluate at the same x nodes"
        )
        if hessian_fit:
            assert ref.replica is None, "Hessian and non-hessian exportgrids mixed"

    # Read the EKO and the operator and theory cards
    eko_op = eko.EKO.read(eko_path)
    theory = eko_op.theory_card
    op = eko_op.operator_card

    _logger.debug(f"Theory card: {json.dumps(theory.raw)}")
    _logger.debug(f"Operator card: {json.dumps(op.raw)}")

    # Check whether the xgrid of the exportgrids matches the EKO
    # if not, construct a copy in memory of the EKO with the xgrid rotated
    # this copy is constructed by replacing the internal inventory of operators
    eko_original_xgrid = eko_op.xgrid
    x_grid = ref.xgrid
    if XGrid(x_grid) != eko_original_xgrid:
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

    all_replicas = []
    for exportgrid in exportgrids:
        pid_to_idx = {pid: idx for idx, pid in enumerate(exportgrid.pids)}
        reindexer = [pid_to_idx[pid] for pid in basis_rotation.flavor_basis_pids]

        all_replicas.append(exportgrid.pdfvalues[reindexer])

    # output is a dictionary {(Q2, nf): (replica, flavour, x)}
    all_evolved, _ = apply.apply_grids(eko_op, np.array(all_replicas))
    # sort the output in terms of (Q2, nf) but grouped by nf
    sorted_evolved = dict(sorted(all_evolved.items(), key=lambda item: (item[0][1], item[0][0])))

    info = info_file.build(theory, op, 1, info_update={})
    info["NumMembers"] = "REPLACE_NREP"
    if hessian_fit:
        info["ErrorType"] = "hessian"
    else:
        info["ErrorType"] = "replicas"
    info["XMin"] = float(x_grid[0])
    info["XMax"] = float(x_grid[-1])
    info["Flavors"] = basis_rotation.flavor_basis_pids
    info.setdefault("NumFlavors", 5)

    return info, sorted_evolved


def evolve_exportgrids_into_lhapdf(eko_path, exportgrids, output_files, info_file):
    """
    Exportgrid evolution function.

    This function takes the path to an ``eko.tar`` and a list of ``ExportGrid`` objects
    and generate the corresponding ``.dat`` LHAPDF files as given by the ``output_files`` list.

    Parameters:
        eko_path: pathlib.Path
            Path to the evolution eko
        exportgrids: list[ExportGrid]
            list of the PDF grids to evolve, the settings must match those of the EKO
        output_files: list[pathlib.Path]
            list of .dat files where the evolved fit will be written, e.g., ["replica_0.dat", "replica_1.dat"]
            the list must be of the same size as exportgrids
        info_file: pathlib.Path
            path to the info file
    """
    if len(exportgrids) != len(output_files):
        raise ValueError("The length of output_files and exportgrids must be equal")

    # Make sure the parent directory of every output_file either exists or can be created
    for f in output_files:
        pathlib.Path(f).parent.mkdir(exist_ok=True, parents=True)
    pathlib.Path(info_file).parent.mkdir(exist_ok=True, parents=True)

    # all evolved is a dictionary {(Q2, nf): (replica, flavour, x)}
    # ordered first by nf and then by Q2 in ascending order
    info, all_evolved = evolve_exportgrid(eko_path, exportgrids)

    # The ``dump`` functions from eko's genpdf are very opinionated regarding the output folder of the files
    # therefore we create a temporary directory where to put stuff and then move it to the right place
    temp_dir = tempfile.TemporaryDirectory()
    temp_path = pathlib.Path(temp_dir.name)

    genpdf.export.dump_info(temp_path, info)
    temp_info = temp_path / f"{temp_path.stem}.info"
    shutil.move(temp_info, info_file)

    # Dump LHAPDF files as .dat files in blocks of nf
    targetgrid = exportgrids[0].xgrid.tolist()
    q2block_per_nf = defaultdict(list)
    for q2, nf in all_evolved.keys():
        q2block_per_nf[nf].append(q2)

    for enum, (exportgrid, output_file) in enumerate(zip(exportgrids, output_files)):
        replica_idx = exportgrid.replica
        blocks = []

        for nf, q2grid in q2block_per_nf.items():

            def pdf_xq2(pid, x, Q2):
                x_idx = targetgrid.index(x)
                pid_idx = info["Flavors"].index(pid)
                ret = x * all_evolved[(Q2, nf)][enum][pid_idx][x_idx]
                return ret

            block = genpdf.generate_block(
                pdf_xq2, xgrid=targetgrid, sorted_q2grid=q2grid, pids=info["Flavors"]
            )
            blocks.append(block)

        dat_path = dump_evolved_replica(blocks, temp_path, replica_idx, exportgrid.hessian)
        shutil.move(dat_path, output_file)

    temp_dir.cleanup()


def evolve_fit(fit_folder, force, eko_path, hessian_fit=False):
    """
    Evolves all the fitted replica in fit_folder/nnfit

    Parameters
    ----------

        fit_folder: str or pathlib.Path
            path to the folder containing the fit
        force: bool
            whether to force the evolution to be done again
        eko_path: str or pathlib.Path
            path where the eko is stored (if None the eko will be
            recomputed)
        hessian_fit: bool
            wether the fit is hessian
    """
    fit_folder = pathlib.Path(fit_folder)
    log_file = fit_folder / LOG_FILE
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

    output_files = []
    exportgrids = []

    for exportgrid_file in fit_folder.glob(f"nnfit/replica_*/{fit_folder.name}.exportgrid"):
        data = yaml_safe.load(exportgrid_file.read_text(encoding="UTF-8"))

        # If hessian fit, force replica to be None
        if hessian_fit:
            data.replica = None

        exportgrids.append(ExportGrid(**data))
        output_files.append(exportgrid_file.with_suffix(".dat"))

    info_path = fit_folder / "nnfit" / f"{fit_folder.name}.info"
    evolve_exportgrids_into_lhapdf(eko_path, exportgrids, output_files, info_path)


def dump_evolved_replica(evolved_blocks, dump_folder, replica_num, hessian_fit=False):
    """
    Dump the evolved replica given by evolved_block as
        dump_folder / f"{dump_folder.stem}_{replica_num:04d}.dat"

    Parameters
    ----------
        evolved_block: list(numpy.array)
            list of blocks of an evolved PDF
        usr_path: pathlib.Path
            path of the fit folder
        replica_num: int
            replica number
        hessian_fit: bool
            wether the fit is hessian
    """
    # create folder to dump the evolved replica if it does not exist
    dump_folder.mkdir(exist_ok=True, parents=True)
    if hessian_fit:
        to_write_in_head = f"PdfType: error\n"
    else:
        to_write_in_head = f"PdfType: replica\nFromMCReplica: {replica_num}\n"
    genpdf.export.dump_blocks(dump_folder, replica_num, evolved_blocks, pdf_type=to_write_in_head)
    return dump_folder / f"{dump_folder.stem}_{replica_num:04d}.dat"
