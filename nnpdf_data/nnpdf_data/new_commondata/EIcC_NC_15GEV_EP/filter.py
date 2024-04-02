from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import yaml

np.random.seed(1234567890)


def read_txt_data(path_txt: Path) -> pd.DataFrame:
    """Return a Panda table in which the columns are ordered as
    [x, Q2, abs]

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


def read_cvs() -> np.ndarray:
    cv_preds = Path("./EIcC_NC_15GEV_EP.yaml")
    cv_yaml = yaml.safe_load(cv_preds.read_text())
    return np.array(cv_yaml["predictions_central"])


def fluctuate_data(central: np.ndarray, abserr: np.ndarray) -> np.ndarray:
    """Fluctuate the central values according to the uncertainties.

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
    df: pd.DataFrame, abserr: Optional[Union[np.ndarray, None]] = None, add_fluctuate: bool = False
) -> None:
    """Write the input kinematics, central values, and uncertainties
    into the new commondata format.

    Parameters
    ----------
    df : pd.DataFrame
        a pandas table containing the information required to generate
        the commondata
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

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # -----------------------------------------------------------------
    # Dump the kinematics
    kins = [
        {
            "x": {"min": None, "mid": float(d["x"]), "max": None},
            "Q2": {"min": None, "mid": float(d["Q2"]), "max": None},
        }
        for _, d in df.iterrows()
    ]

    kinematics_yaml = {"bins": kins}
    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # -----------------------------------------------------------------
    # Dump the uncertainty values
    errors = [
        {"stat": float(d["abs"]), "sys": 0.0, "norm": 2 * c * 1e-2}
        for _, (c, d) in zip(data_central, df.iterrows())
    ]

    error_definition = {
        "stat": {"description": "statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "sys": {"description": "systematic uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "norm": {
            "description": "artificial normalization uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": errors}
    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    input_txt = Path("./EIcC_3_20_A1c_100fb-1.txt")
    df = read_txt_data(input_txt)
    cv_preds = read_cvs()
    fluctuated_cv = fluctuate_data(cv_preds, df["abs"].values)
    write_data(df, abserr=fluctuated_cv, add_fluctuate=True)
