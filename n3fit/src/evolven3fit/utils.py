from scipy.interpolate import interp1d
import numpy as np
import math
import pathlib
import yaml
import shutil


class LhapdfLike:
    pids_dict = {
        -6: "TBAR",
        -5: "BBAR",
        -4: "CBAR",
        -3: "SBAR",
        -2: "UBAR",
        -1: "DBAR",
        21: "GLUON",
        1: "D",
        2: "U",
        3: "S",
        4: "C",
        5: "B",
        6: "T",
        22: "PHT",
    }
    pids_order = [
        "TBAR",
        "BBAR",
        "CBAR",
        "SBAR",
        "UBAR",
        "DBAR",
        "GLUON",
        "D",
        "U",
        "S",
        "C",
        "B",
        "T",
        "PHT",
    ]

    def __init__(self, pdf_grid, q20, x_grid):
        self.pdf_grid = pdf_grid
        self.q20 = q20
        self.x_grid = x_grid
        self.funcs = [
            interp1d(self.x_grid, self.pdf_grid[pid])
            for pid in range(len(self.pids_order))
        ]

    def xfxQ2(self, pid, x, q2):
        if not math.isclose(q2, self.q20, rel_tol=1e-6):
            raise ValueError("The q2 requested is not the fitting scale of this pdf")
        return x * (self.funcs[self.pids_order.index(self.pids_dict[pid])](x))

    def hasFlavor(self, pid):
        if pid in self.pids_dict.keys():
            return True
        return False


def read_runcard(usr_path):
    """
    reads the runcard and returns the relevant information for evolven3fit
    """
    runcard_path = usr_path / "filter.yml"
    with runcard_path.open() as fp:
        data = yaml.safe_load(fp)
    return data


def generate_q2grid(Q0, Qfin):
    """
    Generate the q2grid used in the final evolved pdfs (Temporary solution)
    """
    grid = np.geomspace(Q0 ** 2, Qfin ** 2, num=100)
    return grid


def mv_file(file_path, dest_path):
    shutil.move(str(file_path), str(dest_path))


def fix_info_path(usr_path):
    info_file_path = usr_path / "nnfit" / usr_path.stem / (usr_path.stem + ".info")
    dest_path_info = usr_path / "nnfit" / (usr_path.stem + ".info")
    mv_file(info_file_path, dest_path_info)


def fix_replica_path(usr_path, replica_num):
    replica_file_path = (
        usr_path / "nnfit" / usr_path.stem / f"{usr_path.stem}_{replica_num:04d}.dat"
    )
    dest_path_replica = (
        usr_path / "nnfit" / f"replica_{replica_num}" / f"{usr_path.stem}.dat"
    )
    mv_file(replica_file_path, dest_path_replica)
