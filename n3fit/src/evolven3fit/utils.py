import pathlib

import numpy as np
import warnings
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


def Q2_to_t(q2: float, lambda2) -> float:
    """Map to define loglog spacing of the Q2grid"""
    return np.log(np.log(q2 / lambda2))

def t_to_Q2(t: float, lambda2) -> float:
    """Map to go from the loglog spacing back to Q2"""
    return lambda2 * np.exp(np.exp(t))

def generate_q2grid(Q0, Qmin, Qmax, match_dict, total_points, total_points_ic, legacy40=False):
    """Generate the q2grid used in the final evolved pdfs or use the default grid if legacy40 is set. 
    The grid uses log(log(Q2/Lambda2)) spacing between the points, with Lambda2=0.0625.

    The grid from Q2_min --> mc is made separately from the rest to be able to let it contain at
    least 5 points. The rest of the grids contains the same loglog spacing from threshold to threshold.

    Q2_min --> mc^2 --> Q0^2 --> mb^2 --> Q2_max

    match_dict contains the quark mass thresholds and factors. Factor is the number to be 
    multiplied to mass in order to obtain the relative matching scale.
    """
    
    # If flag --legacy40 is set return handmade legacy grid
    if legacy40:
        return Q2GRID_NNPDF40
    # Otherwise dynamically create the grid from Q2_min --> Q2_max
        Q2_min = Qmin**2 # 1.0**2
        Q2_max = Qmax**2 # 1e5**2                    
        LAMBDA2 = 0.0625                  
        
        # Collect all node Q2's from Q0^2 --> Q2_max
        q0_2 = Q0**2
        node_Q2 = [q0_2, (match_dict["mb"]*match_dict["kbThr"])**2, Q2_max]
            
        # Make initial uniform grid in t from Q0^2 --> Q2_max
        t_min = Q2_to_t(q0_2, LAMBDA2)
        t_max = Q2_to_t(Q2_max, LAMBDA2)
        t_vals = np.linspace(t_min, t_max, total_points)
        q2_vals = t_to_Q2(t_vals, LAMBDA2)
        
        # Count how many points fall into each subgrid
        n_intervals = len(node_Q2) - 1
        nQpoints = np.zeros(n_intervals, dtype=int)
        subgridindex = 0
            
        for q2 in q2_vals:
            while subgridindex < n_intervals - 1 and q2 >= node_Q2[subgridindex + 1]:
                subgridindex += 1
            nQpoints[subgridindex] += 1
        
        # Make t grid from Q2_min --> Q0^2
        t_min_ic = Q2_to_t(Q2_min, LAMBDA2)
        t_max_ic = Q2_to_t((match_dict["mc"]*match_dict["kcThr"])**2, LAMBDA2)
        t_vals_ic = np.linspace(t_min_ic, t_max_ic, total_points_ic)
        q2_vals_ic = t_to_Q2(t_vals_ic, LAMBDA2)

        # Now build each subgrid to contain the points we want
        grids = []
        grids.append(q2_vals_ic)
        for i in range(len(node_Q2) - 1):
            q2_lo, q2_hi = node_Q2[i], node_Q2[i + 1]
            t_lo, t_hi = Q2_to_t(q2_lo, LAMBDA2), Q2_to_t(q2_hi, LAMBDA2)
            npts = int(nQpoints[i])
            t_subgr = np.linspace(t_lo, t_hi, npts)
            q2_subgr = t_to_Q2(t_subgr, LAMBDA2)
            if i < n_intervals - 1:
                q2_subgr = q2_subgr[:-1]
            grids.append(q2_subgr)

        # Combine all subgrids and return as an array
        q2_full = np.concatenate(grids)
        q2_full[0] = Q2_min
        q2_full[-1] = Q2_max
        return q2_full
    
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
