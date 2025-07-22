import pathlib

import numpy as np

from validphys.utils import yaml_safe

from .q2grids import Q2GRID_DEFAULT, Q2GRID_NNPDF40


def read_runcard(usr_path):
    """Read the runcard and return the relevant information for evolven3fit"""
    return yaml_safe.load((usr_path / "filter.yml").read_text(encoding="UTF-8"))


def get_theoryID_from_runcard(usr_path):
    """Return the theoryID from the runcard"""
    # read the runcard
    my_runcard = read_runcard(usr_path)
    return my_runcard["theory"]["theoryid"]


def generate_q2grid(Q0, Qfin, Q_points, match_dict, nf0=None, legacy40=False, theory_41=False):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if Qfin or Q_points is
    not provided.

    match_dict contains the couples (mass : factor) where factor is the number to be multiplied to mass
    in order to obtain the relative matching scale.
    """
    if Qfin is None and Q_points is None:
        if legacy40:
            return Q2GRID_NNPDF40
        elif nf0 in (3, 4, 5):
            return Q2GRID_DEFAULT
        elif theory_41:
            return Q2GRID_DEFAULT_41
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
                grids.append(np.geomspace(Q_ini**2, match_scale**2, num=num_points, endpoint=False))
                Q_ini = match_scale
        num_points = Q_points - sum(num_points_list)
        grids.append(np.geomspace(Q_ini**2, Qfin**2, num=num_points))
        return np.concatenate(grids).tolist()


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


def check_filter(config_folder):
    """Check if config_folder contains a filter.yml file."""
    usr_path = pathlib.Path(config_folder)
    filter_path = usr_path / "filter.yml"
    if not filter_path.is_file():
        raise ValueError(
            f"filter.yaml file not found: the path {filter_path.absolute()} is not valid"
        )


def check_nnfit_folder(config_folder):
    """Check if config_folder contains a nnfit folder."""
    usr_path = pathlib.Path(config_folder)
    nnfitpath = usr_path / "nnfit"
    if not nnfitpath.is_dir():
        raise ValueError("nnfit folder not found: provided path is not valid")
