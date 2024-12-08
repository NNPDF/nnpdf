from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

# Definition of various topologies used in Polarized Dijets
# NOTE: the observable is symmetric for jet1 and jet2,
# so 1 and 2 are not ordered in pT.
TOPO_DEF = {
    "A": {"abs_eta1_min": 0.3, "abs_eta1_max": 0.9, "abs_eta2_min": 0.3, "abs_eta2_max": 0.9},
    "B": {"abs_eta1_min": 0, "abs_eta1_max": 0.3, "abs_eta2_min": 0.3, "abs_eta2_max": 0.9},
    "C": {"abs_eta1_min": 0, "abs_eta1_max": 0.3, "abs_eta2_min": 0, "abs_eta2_max": 0.3},
    "D": {"abs_eta1_min": 0.3, "abs_eta1_max": 0.9, "abs_eta2_min": 0.3, "abs_eta2_max": 0.9},
    "I": {"abs_eta_min": 0, "abs_eta_max": 0.9},
}


def read_central_values(path: Path) -> np.ndarray:
    """Read the central values from the theory predictions.

    Parameters
    ----------
    path: Path
        path to the file containing the central values

    Returns
    -------
    np.ndarray:
        array containing the central predictions
    """
    cv_yaml = yaml.safe_load(path.read_text())
    return np.array(cv_yaml["predictions_central"])


def read_txt_data(path_txt: Path) -> pd.DataFrame:
    """Reads a rawdata `.txt` file and returns a pandas dataframe in which the
    columns are ordered as [x, Q2, abs]

    Parameters
    ----------
    path_txt : Path
        path to the .txt file

    Returns
    -------
    pd.DataFrame
        a panda table containing the values
    """
    colnames = ["x", "Q2", "abs"]
    return pd.read_csv(
        path_txt, delim_whitespace=True, names=colnames, usecols=[i for i in range(len(colnames))]
    )


def read_excel(path_xlsx: Path, beams: tuple) -> pd.DataFrame:
    """Parse the xlsx file containing all the information regarding
    the projections and returns the ones corresponding to the chosen
    beam energies.

    Parameters
    ----------
    path_xlsx : Path
        path to the xlsx file
    beams: tuple
        tuple specifying the beam energies of the lepton & proton

    Returns
    -------
    pd.DataFrame
        returns a panda table corresponding to the chosen beams
    """
    xdf = pd.read_excel(path_xlsx)
    el, ep = beams
    df_beam = xdf[(xdf["El"] == el) & (xdf["Eh"] == ep)]
    # Because sometimes the raw data contains events with `y>1` due to the
    # detector and MC simulation, we need to make sure that these are not
    # considered and removed.
    return df_beam[df_beam["y"] <= 1]


def fluctuate_data(central: np.ndarray, abserr: np.ndarray) -> np.ndarray:
    """Fluctuate the central values according to the uncertainties. This
    is done because when the measurements are not real data but rather
    pseudodata projections, the central values are the ones computed from
    theory predictions and fluctuaded according to the projected uncertainties.

    Parameters
    ----------
    central : np.ndarray
        array of central values
    abserr : np.ndarray
        array containing the values of the errors

    Returns
    -------
    np.ndarray
        fluctuated central values according to a normal distribution
    """
    shifted_cv = [np.random.normal(c, e) for c, e in zip(central, abserr)]
    return np.array(shifted_cv)


def write_data(
    df: pd.DataFrame,
    asym: bool = True,
    abserr: Optional[Union[np.ndarray, None]] = None,
    add_fluctuate: bool = False,
) -> None:
    """Write the input kinematics, central values, and uncertainties
    into the new commondata format.

    Parameters
    ----------
    df : pd.DataFrame
        a pandas table containing the information required to generate
        the commondata
    asym: bool
        boolean that states whether or not the observable is an Asymmetry,
        as is the case for ATHENA
    abserr: Optional[Union[np.ndarray, None]]
        if not None contains the fluctuated centra values
    add_fluctuate: bool
        whether or not to fluctuate the central values
    """
    # -----------------------------------------------------------------
    # Dump the Central values
    if not add_fluctuate:
        data_central = [None for _ in range(len(df))]
    else:
        data_central = abserr.tolist()
    print(f"The dataset has {len(data_central)} datapoints!")

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    if asym:
        # -----------------------------------------------------------------
        # Prepare the kinematics
        kins = [
            {
                "x": {"min": None, "mid": float(d["x"]), "max": None},
                "Q2": {"min": None, "mid": float(d["Q2"]), "max": None},
                "y": {"min": None, "mid": float(d["y"]), "max": None},
            }
            for _, d in df.iterrows()
        ]

        # -----------------------------------------------------------------
        # Prepare the uncertainty values
        errors = []
        for idx, (_, d) in enumerate(df.iterrows()):
            if not add_fluctuate:
                errors.append({"stat": None, "sys": None, "shift_lumi": None, "norm": None})
            else:
                errors.append(
                    {
                        "stat": float(data_central[idx] * d["unpol_stat_percent"] * 1e-2),
                        "sys": float(data_central[idx] * d["ptpt_percent"] * 1e-2),
                        "shift_lumi": float(d["shift_uncer"]),
                        "norm": float(data_central[idx] * d["norm_percent"] * 1e-2),
                    }
                )

        error_definition = {
            "stat": {
                "description": "statistical uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "sys": {"description": "systematic uncertainty", "treatment": "MULT", "type": "UNCORR"},
            "shift_lumi": {
                "description": "uncertainty on the precision of the relative luminosity",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "norm": {
                "description": "relative (percent) normalization uncertainty (beam pol)",
                "treatment": "MULT",
                "type": "CORR",
            },
        }
    else:
        # -----------------------------------------------------------------
        # Prepare the kinematics
        kins = [
            {
                "x": {"min": None, "mid": float(d["x"]), "max": None},
                "Q2": {"min": None, "mid": float(d["Q2"]), "max": None},
            }
            for _, d in df.iterrows()
        ]

        # -----------------------------------------------------------------
        # Prepare the uncertainty values
        errors = [{"stat": float(d["abs"]), "sys": 0.0} for _, d in df.iterrows()]

        error_definition = {
            "stat": {
                "description": "statistical uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "sys": {"description": "systematic uncertainty", "treatment": "ADD", "type": "UNCORR"},
        }

    # -----------------------------------------------------------------
    # Dump the kinematics and uncertainties
    kinematics_yaml = {"bins": kins}
    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    uncertainties_yaml = {"definitions": error_definition, "bins": errors}
    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)
