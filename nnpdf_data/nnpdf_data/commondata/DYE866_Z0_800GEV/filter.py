from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import typing
from typing import List

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.hera_utils import commondata, covmat_is_close


def readdata() -> pd.DataFrame:
    col_names = [
        "xF",
        "Mmin",
        "Mmax",
        "Mavg",
        "xFavg",
        "pt",
        "dsig",
        "stat",
        "syst",
        "kfact",
        "rsig",
        "rstat",
        "rsyst",
    ]
    table_path = Path(f"./rawdata/table.csv")
    df = pd.read_csv(table_path, names=col_names)
    return df


@dataclass
class E866commondata(commondata):
    def __init__(self, data: pd.DataFrame, dataset_name: str, process: str):

        # Definitions, compute Jacobian, get dsig/dy/dM
        M = (data["Mmax"] + data["Mmin"]) / 2
        M2 = M * M
        sqrts = M / M * 38.8
        s = sqrts**2
        tau = M**2 / s
        tau = tau.to_numpy()
        xF = data["xF"]
        y = np.arcsinh(xF / np.sqrt(tau) / 2)
        jac = np.sqrt(xF**2 + 4 * tau)
        dsigdydM = data["dsig"] * jac

        # Set the central values
        self.central_values = dsigdydM.astype(float).to_numpy()

        # Pick the the kinematic quantities
        kin = pd.concat([y, M2, sqrts], axis=1)
        kin = kin.set_axis(["y", "M2", "sqrts"], axis=1)
        self.kinematics = kin.astype(float).to_numpy()
        self.kinematic_quantities = ["y", "M2", "sqrts"]

        # Statistical uncertainties.
        self.statistical_uncertainties = data["stat"] * jac

        # Systematic uncertainty
        syst = data["syst"] * jac

        # Normalisation uncertainty of 6.5% from beam intensity calibration.
        norm = 6.5 / 100
        norm = norm * self.central_values

        self.systematic_uncertainties = np.dstack((syst, norm))[0]
        self.systypes = [("ADD", "UNCORR"), ("MULT", "CORR")]

        self.process = process
        self.dataset_name = dataset_name


def main():
    data = readdata()
    # First create the commondata variant without the nuclear uncertainties.
    DYE866 = E866commondata(data, "DYE866_Z0", "Z0")
    DYE866.write_new_commondata(
        Path("data_reimplemented_PXSEC.yaml"),
        Path("kinematics_reimplemented_PXSEC.yaml"),
        Path("uncertainties_reimplemented_PXSEC.yaml"),
    )

    if covmat_is_close("DYE866_Z0_800GEV_PXSEC", "legacy", "reimplemented"):
        print("covmat is close")
    else:
        print("covmat is different.")


if __name__ == "__main__":
    main()
