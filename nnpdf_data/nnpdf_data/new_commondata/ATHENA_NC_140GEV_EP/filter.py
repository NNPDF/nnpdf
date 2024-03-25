import yaml
import numpy as np
import pandas as pd

from pathlib import Path
from typing import Optional, Union

np.random.seed(1234567890)


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
    return xdf[(xdf["El"] == el) & (xdf["Eh"] == ep)]


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
    df: pd.DataFrame,
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
        raise ValueError("This is not supported yet!")

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)
    print(f"[+] The number of data points is {len(data_central)}")

    # -----------------------------------------------------------------
    # Dump the kinematics
    kins = [
        {
            "x": {"min": None, "mid": float(d["x"]), "max": None},
            "Q2": {"min": None, "mid": float(d["Q2"]), "max": None},
            "y": {"min": None, "mid": float(d["y"]), "max": None},
        }
        for _, d in df.iterrows()
    ]

    kinematics_yaml = {"bins": kins}
    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)
    print(f"[+] The number of kinematic points is {len(kins)}")

    # -----------------------------------------------------------------
    # Dump the uncertainty values
    errors = []
    for idx, (_, d) in enumerate(df.iterrows()):
        if not add_fluctuate:
            errors.append(
                {"stat": None, "sys": None, "shift_lumi": None, "norm": None}
            )
        else:
            errors.append(
                {
                    "stat": data_central[idx] * df["unpol_stat_percent"] * 1e-2,
                    "sys": data_central[idx] * df["ptpt_percent"] * 1e-2,
                    "shift_lumi": df["shift_uncer"],
                    "norm": data_central[idx] * df["norm_percent"] * 1e-2,
                }
            )
    print(f"[+] The number of uncertainty points is {len(errors)}")

    error_definition = {
        "stat": {
            "description": "statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "sys": {
            "description": "systematic uncertainty",
            "treatment": "MULT", # TODO: to check
            "type": "UNCORR",
        },
        "shift_lumi": {
            "description": "uncertainty on the precision of the relative luminosity",
            "treatment": "ADD",
            "type": "UNCORR", # TODO: to check
        },
        "norm": {
            "description": "relative (percent) normalization uncertainty (beam pol)",
            "treatment": "MULT", # TODO: to check
            "type": "CORR", # TODO: to check
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": errors}
    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    BEAMS = (18, 275)
    input_xlsx = Path("./ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    write_data(xdf)
