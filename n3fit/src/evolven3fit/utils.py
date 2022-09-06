from scipy.interpolate import interp1d
import numpy as np
import math
import pathlib
from reportengine.compat import yaml
import shutil


class LhapdfLike:
    """
    Class which emulates lhapdf but only for an initial condition PDF (i.e. with only one q2 value)
    """

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
            interp1d(self.x_grid, self.pdf_grid[pid], kind="cubic")
            for pid in range(len(self.pids_order))
        ]

    def xfxQ2(self, pid, x, q2):
        """
        Return the value of the PDF for the requested pid and x value. If the requested q2
        value is different from the (only) value available, it raises an error.

        Parameters
        ----------

            pid: int
                pid index of particle
            x: float
                x-value
            q2: float
                Q square value

        Returns
        -------
            : float
            x * PDF value
        """
        if not math.isclose(q2, self.q20, rel_tol=1e-6):
            raise ValueError("The q2 requested is not the fitting scale of this pdf")
        return self.funcs[self.pids_order.index(self.pids_dict[pid])](x)

    def hasFlavor(self, pid):
        """
        Check if the requested pid is in the PDF.
        """
        return pid in self.pids_dict.keys()


def read_runcard(usr_path):
    """
    reads the runcard and returns the relevant information for evolven3fit
    """
    return yaml.safe_load((usr_path / "filter.yml"))


def generate_q2grid(Q0, Qfin):
    """
    Generate the q2grid used in the final evolved pdfs (Temporary solution)
    """
    return np.geomspace(Q0**2, Qfin**2, num=20).tolist()


def generate_x_grid():
    """
    Generate the xgrid used for the eko
    """
    grid = np.geomspace(1e-09, 1.0, num=196).tolist()
    return grid

def fix_info_path(usr_path):
    """Fix the location of the info file from the folder nnfit/usr_path to 
    just nnfit
    """
    nnfit = usr_path / "nnfit"
    info_file = usr_path.stem + ".info"
    info_file_path = nnfit / usr_path.stem / info_file
    dest_path_info = nnfit / info_file
    shutil.move(info_file_path, dest_path_info)


def fix_replica_path(usr_path, replica_num):
    """Fix the location of the dat file of the replica <replica_num> from the folder nnfit/usr_path to 
    just nnfit/replica_<replica_num>
    """
    nnfit = usr_path / "nnfit"
    replica_file_path = nnfit / usr_path.stem / f"{usr_path.stem}_{replica_num:04d}.dat"
    dest_path_replica = nnfit / f"replica_{replica_num}" / f"{usr_path.stem}.dat"
    shutil.move(replica_file_path, dest_path_replica)
