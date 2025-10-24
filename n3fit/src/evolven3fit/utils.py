import pathlib

import numpy as np

from validphys.utils import yaml_safe

from .q2grids import Q2GRID_NNPDF40


def read_runcard(usr_path):
    """Read the runcard and return the relevant information for evolven3fit"""
    return yaml_safe.load((usr_path / "filter.yml").read_text(encoding="UTF-8"))


def get_theoryID_from_runcard(usr_path):
    """Return the theoryID from the runcard"""
    # read the runcard
    my_runcard = read_runcard(usr_path)
    return my_runcard["theory"]["theoryid"]


def generate_q2grid(Q0, Qfin, Q_points, match_dict, legacy40=False):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if Qfin or Q_points is
    not provided.

    match_dict contains the couples (mass : factor) where factor is the number to be multiplied to mass
    in order to obtain the relative matching scale.
    """
    if Qfin is None and Q_points is None:
        if legacy40:
            return Q2GRID_NNPDF40
        else: 
            Q2_min = 1.0**2
            Q2_max = 1e5**2
            Q0 = 1.65                    
            Lambda2 = 0.0625                  
            total_points = 50
            min_per_batch = 2

            # Collect threshold Q2's
            node_Q2 = []
            q0_2 = Q0**2
            node_Q2.append(q0_2)
            node_Q2.append(match_dict["mb"]**2)
            node_Q2.append(Q2_max)
            node_Q2 = sorted(set(node_Q2))

            # Define map
            def Q2_to_t(q2: float) -> float:
                return np.log(np.log(q2 / Lambda2))

            def t_to_Q2(t: float) -> float:
                return Lambda2 * np.exp(np.exp(t))
            
            # Make initial uniform grid in t from Q0^2 to Q2_max
            t_min = Q2_to_t(q0_2)
            t_max = Q2_to_t(Q2_max)
            t_vals = np.linspace(t_min, t_max, total_points)
            q2_vals = t_to_Q2(t_vals)

            # Make t grid from Q2_min t Q0^2
            t_vals_ic = np.linspace(Q2_to_t(Q2_min), Q2_to_t(match_dict["mc"]**2), 6)
            q2_vals_ic = t_to_Q2(t_vals_ic)
            
            # Count how many points fall into each subgrid
            n_intervals = len(node_Q2) - 1
            nQpoints = np.zeros(n_intervals, dtype=int)
            subgridindex = 0
            
            for q2 in q2_vals:
                while subgridindex < n_intervals - 1 and q2 >= node_Q2[subgridindex + 1]:
                    subgridindex += 1
                nQpoints[subgridindex] += 1

            # Now build each subgrid to contain the points we want
            grids = []
            grids.append(q2_vals_ic)
            for i in range(n_intervals):
                q2_lo, q2_hi = node_Q2[i], node_Q2[i + 1]
                t_lo, t_hi = Q2_to_t(q2_lo), Q2_to_t(q2_hi)
                npts = int(nQpoints[i])
                t_subgr = np.linspace(t_lo, t_hi, npts)
                q2_subgr = t_to_Q2(t_subgr)
                if i < n_intervals - 1:
                    q2_subgr = q2_subgr[:-1]
                grids.append(q2_subgr)
            
            # Combine all subgrids
            q2_full = np.concatenate(grids)
            q2_full[0] = Q2_min
            q2_full[-1] = Q2_max
            return np.sqrt(q2_full.tolist())      
            
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
