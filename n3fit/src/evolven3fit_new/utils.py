import pathlib
import shutil

import numpy as np
from scipy.interpolate import interp1d

from reportengine.compat import yaml
from validphys.pdfbases import PIDS_DICT

Q2GRID_Nf04 = (
    np.array(
        [
            1.6500000e00,
            1.7874388e00,
            1.9429053e00,
            2.1193749e00,
            2.3204100e00,
            2.5502944e00,
            2.8142025e00,
            3.1184122e00,
            3.4705775e00,
            3.8800751e00,
            4.3584516e00,
            4.9200000e00,
            5.5493622e00,
            6.2897452e00,
            7.1650687e00,
            8.2052867e00,
            9.4481248e00,
            1.0941378e01,
            1.2745972e01,
            1.4940062e01,
            1.7624572e01,
            2.0930715e01,
            2.5030298e01,
            3.0149928e01,
            3.6590777e01,
            4.4756282e01,
            5.5191298e01,
            6.8637940e01,
            8.6115921e01,
            1.0903923e02,
            1.3938725e02,
            1.7995815e02,
            2.3474820e02,
            3.0952544e02,
            4.1270732e02,
            5.5671861e02,
            7.6011795e02,
            1.0509694e03,
            1.4722574e03,
            2.0906996e03,
            3.0112909e03,
            4.4016501e03,
            6.5333918e03,
            9.8535186e03,
            1.5109614e04,
            2.3573066e04,
            3.7444017e04,
            6.0599320e04,
            1.0000000e05,
        ]
    )
    ** 2
)

Q2GRID_Nf03 = (
    np.array(
        [
            1.0000000e00,
            1.0768843e00,
            1.1642787e00,
            1.2640247e00,
            1.3783565e00,
            1.5100000e00,
            1.6573843e00,
            1.8279487e00,
            2.0263188e00,
            2.2582323e00,
            2.5308507e00,
            2.8531703e00,
            3.2365690e00,
            3.6955380e00,
            4.2486693e00,
            4.9200000e00,
            5.6571821e00,
            6.5475141e00,
            7.6300446e00,
            8.9555329e00,
            1.0590474e01,
            1.2622686e01,
            1.5169120e01,
            1.8386905e01,
            2.2489085e01,
            2.7767274e01,
            3.4624624e01,
            4.3624282e01,
            5.5561424e01,
            7.1571582e01,
            9.3295496e01,
            1.2313315e02,
            1.6464038e02,
            2.2315640e02,
            3.0681103e02,
            4.2816505e02,
            6.0692308e02,
            8.7449251e02,
            1.2817733e03,
            1.9127020e03,
            2.9082314e03,
            4.5095982e03,
            7.1379509e03,
            1.1543948e04,
            1.9094934e04,
            3.2338760e04,
            5.6137084e04,
            1.0000000e05,
        ]
    )
    ** 2
)


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
    """Read the runcard and return the relevant information for evolven3fit_new"""
    return yaml.safe_load((usr_path / "filter.yml").read_text(encoding="UTF-8"))


def get_theoryID_from_runcard(usr_path):
    """Return the theoryID from the runcard"""
    # read the runcard
    my_runcard = read_runcard(usr_path)
    return my_runcard["theory"]["theoryid"]


def generate_q2grid(Q0, Qfin, Q_points, match_dict, nf0=None):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if Qfin or Q_points is
    not provided.

    match_dict contains the couples (mass : factor) where factor is the number to be multiplied to mass
    in order to obtain the relative matching scale.
    """
    if Qfin is None and Q_points is None:
        if nf0 == 4:
            return Q2GRID_Nf04
        elif nf0 == 3:
            return Q2GRID_Nf03
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
