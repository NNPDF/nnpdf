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
        """Return the value of the PDF for the requested pid and x value. If the requested q2
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
        """Check if the requested pid is in the PDF.
        """
        return pid in self.pids_dict.keys()


def read_runcard(usr_path):
    """Read the runcard and return the relevant information for evolven3fit
    """
    return yaml.safe_load((usr_path / "filter.yml").read_text())


def generate_q2grid(Q0, Qfin, Q_points, match_dict):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if Qfin or Q_points is
    not provided.
    """
    if (Qfin is None) or (Q_points is None):   
        return (np.array([1.6500000e+00, 1.7874388e+00, 1.9429053e+00, 2.1193749e+00, 2.3204100e+00,
                2.5502944e+00, 2.8142025e+00, 3.1184122e+00, 3.4705775e+00, 3.8800751e+00, 4.3584516e+00,
                4.9200000e+00, 4.9200000e+00, 5.5493622e+00, 6.2897452e+00, 7.1650687e+00, 8.2052867e+00,
                9.4481248e+00, 1.0941378e+01, 1.2745972e+01, 1.4940062e+01, 1.7624572e+01, 2.0930715e+01,
                2.5030298e+01, 3.0149928e+01, 3.6590777e+01, 4.4756282e+01, 5.5191298e+01, 6.8637940e+01,
                8.6115921e+01, 1.0903923e+02, 1.3938725e+02, 1.7995815e+02, 2.3474820e+02, 3.0952544e+02,
                4.1270732e+02, 5.5671861e+02, 7.6011795e+02, 1.0509694e+03, 1.4722574e+03, 2.0906996e+03,
                3.0112909e+03, 4.4016501e+03, 6.5333918e+03, 9.8535186e+03, 1.5109614e+04, 2.3573066e+04,
                3.7444017e+04, 6.0599320e+04, 1.0000000e+05])**2)
    grids = []
    Q_ini=Q0
    num_points_list = []
    for masses in match_dict.keys():
        match_scale = masses*match_dict[masses]
        num_points = int(Q_points*(np.log(match_scale/Q0)/np.log(Qfin/Q_ini)))
        num_points_list.append(num_points)
        grids.append(np.geomspace(Q0**2, match_scale**2, num=num_points))
        Q0 = match_scale
    num_points = Q_points - sum(num_points_list)
    grids.append(np.geomspace(Q0**2, Qfin**2, num=num_points ))
    return np.concatenate(grids).tolist()

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
