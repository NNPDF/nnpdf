import pathlib
import shutil

import numpy as np
from scipy.interpolate import interp1d

from reportengine.compat import yaml
from validphys.pdfbases import PIDS_DICT

from .q2grids import Q2GRID_DEFAULT, Q2GRID_NNPDF40


class LhapdfLike:
    """
    Class which emulates lhapdf but only for an initial condition PDF (i.e. with only one q2 value).

    Q20 is the fitting scale fo the pdf and it is the only available scale for the objects of this class.

    X_GRID is the grid of x values on top of which the pdf is interpolated.

    PDF_GRID is a dictionary containing the pdf grids at fitting scale for each pid.
    """

    def __init__(self, pdf_grid, q20, x_grid):
        self.pdf_grid = pdf_grid
        self.q20 = q20
        self.x_grid = x_grid
        self.funcs = [
            interp1d(self.x_grid, self.pdf_grid[pid], kind="cubic") for pid in range(len(PIDS_DICT))
        ]

    def xfxQ2(self, pid, x, q2):
        """Return the value of the PDF for the requested pid, x value and, whatever the requested
        q2 value, for the fitting q2.

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
        return self.funcs[list(PIDS_DICT.values()).index(PIDS_DICT[pid])](x)

    def hasFlavor(self, pid):
        """Check if the requested pid is in the PDF."""
        return pid in PIDS_DICT


def read_runcard(usr_path):
    """Read the runcard and return the relevant information for evolven3fit"""
    return yaml.safe_load((usr_path / "filter.yml").read_text(encoding="UTF-8"))


def get_theoryID_from_runcard(usr_path):
    """Return the theoryID from the runcard"""
    # read the runcard
    my_runcard = read_runcard(usr_path)
    return my_runcard["theory"]["theoryid"]


def generate_q2grid(Q0, Qfin, Q_points, match_dict, nf0=None, legacy40=False):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if Qfin or Q_points is
    not provided.

    match_dict contains the couples (mass : factor) where factor is the number to be multiplied to mass
    in order to obtain the relative matching scale.
    """
    if Qfin is None and Q_points is None:
        if legacy40:
            return Q2GRID_NNPDF40
        elif nf0 in (3, 4):
            return Q2GRID_DEFAULT
        elif nf0 is None:
            raise ValueError("In order to use a default grid, a value of nf0 must be provided")
        else:
            raise NotImplementedError(f"No default grid in Q available for {nf0=}")
    elif Qfin is None or Q_points is None:
        raise ValueError("q_fin and q_points must be specified either both or none of them")
    else:
        grids = []
        Q_ini = Q0
        num_points_list = []
        for masses in match_dict:
            match_scale = masses * match_dict[masses]
            # Fraction of the total points to be included in this batch is proportional
            # to the log of the ratio between the initial scale and final scale of the
            # batch itself (normalized to the same log of the global initial and final
            # scales)
            if match_scale < Qfin:
                frac_of_point = np.log(match_scale / Q_ini) / np.log(Qfin / Q0)
                num_points = int(Q_points * frac_of_point)
                num_points_list.append(num_points)
                grids.append(
                    np.geomspace(Q_ini**2, match_scale**2, num=num_points, endpoint=False)
                )
                Q_ini = match_scale
        num_points = Q_points - sum(num_points_list)
        grids.append(np.geomspace(Q_ini**2, Qfin**2, num=num_points))
        return np.concatenate(grids).tolist()


def fix_info_path(usr_path):
    """Fix the location of the info file from the folder nnfit/usr_path to
    just nnfit

    Examples
    --------
    Starting from the info path
        initial_info_file_path = "/myfolder/myfit/nnfit/myfit/myfit.info"
    and using this function with usr_path = "/myfolder/myfit", one gets
        final_info_file_path = "/myfolder/myfit/nnfit/myfit.info"
    """
    nnfit = usr_path / "nnfit"
    info_file = usr_path.stem + ".info"
    info_file_path = nnfit / usr_path.stem / info_file
    dest_path_info = nnfit / info_file
    shutil.move(info_file_path, dest_path_info)


def fix_replica_path(usr_path, replica_num):
    """Fix the location of the dat file of the replica <replica_num> from the folder nnfit/usr_path to
    just nnfit/replica_<replica_num>

    Examples
    --------
    Starting from the replica 5 path
        initial_replica_file_path = "/myfolder/myfit/nnfit/myfit/myfit_5.dat"
    and using this function with usr_path = "/myfolder/myfit", one gets
        final_replica_file_path = "/myfolder/myfit/nnfit/replica_5/myfit.dat"
    """
    nnfit = usr_path / "nnfit"
    replica_file_path = nnfit / usr_path.stem / f"{usr_path.stem}_{replica_num:04d}.dat"
    dest_path_replica = nnfit / f"replica_{replica_num}" / f"{usr_path.stem}.dat"
    shutil.move(replica_file_path, dest_path_replica)


def check_is_a_fit(config_folder):
    """Check if config_folder is a fit folder, i.e. if it contains the filter.yml file
    and the nnfit folder."""
    usr_path = pathlib.Path(config_folder)
    filter_path = usr_path / "filter.yml"
    if not filter_path.is_file():
        raise ValueError(
            "filter.yaml file not found: the path" + str(filter_path.absolute()) + " is not valid"
        )
    nnfitpath = usr_path / "nnfit"
    if not nnfitpath.is_dir():
        raise ValueError("nnfit folder not found: provided path is not valid")
