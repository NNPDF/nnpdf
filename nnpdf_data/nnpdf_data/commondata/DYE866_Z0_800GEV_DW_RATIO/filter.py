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
    col_names = ["x2", "xlow", "xhigh", "xf", "pt", "M", "sig", "errp", "errm"]
    table_path = Path(f"./rawdata/table11.csv")
    df = pd.read_csv(table_path, header=12, names=col_names, nrows=15)
    return df


def nuclear_uncert_dw(tableN: PathLike, tablep: PathLike):
    dfN = pd.read_table(tableN)
    dfp = pd.read_table(tablep)
    return dfN, dfp


@dataclass
class E866_DW_RATIO_commondata(commondata):
    def __init__(self, data: pd.DataFrame, dataset_name: str, process: str):

        # Kinematic quantities.
        self.central_values = data["sig"].astype(float).to_numpy()
        x2 = 0.5 * (data["xhigh"] + data["xlow"])
        x1 = data["xf"] + data["x2"]
        s = x2 / x2 * 38.8 * 38.8
        tau = data["M"] ** 2 / 38.8**2
        y = np.log((data["M"] / np.sqrt(s)) / data["x2"])
        kin = pd.concat([y, data["M"] ** 2, np.sqrt(s)], axis=1)
        kin = kin.set_axis(["y", "M2", "s"], axis=1)
        self.kinematics = kin.astype(float).to_numpy()
        self.kinematic_quantities = ["y", "M2", "sqrts"]

        # Statistical uncertainties.
        self.statistical_uncertainties = data["errp"].astype(float).to_numpy()

        # Systematic uncertainties. Taken from Phys. Rev. D. Vol. 64, 052002 table 10
        syst = 0.97 / 100
        syst = pd.DataFrame(syst * data["sig"])
        systypes = [("ADD", "CORR")]

        # Compute the point-to-point uncertainties
        nrep = 100
        norm = np.sqrt(nrep)
        dfN, dfp = nuclear_uncert_dw(
            "rawdata/nuclear_ite/output/tables/group_result_table.csv",
            "rawdata/proton_ite/output/tables/group_result_table.csv",
        )
        for rep in range(1, nrep + 1):
            Delta = (dfN[f"rep_{rep:05d}"] - dfp["theory_central"]) / norm
            syst[f"NUCLEAR{rep:03d}"] = Delta
            systypes.append(("ADD", f"NUCLEAR{rep:03d}"))

        self.systematic_uncertainties = syst.astype(float).to_numpy()
        self.systypes = systypes
        self.process = process
        self.dataset_name = dataset_name


def main():
    data = readdata()
    # First create the commondata variant without the nuclear uncertainties.
    DYE866 = E866_DW_RATIO_commondata(data, "DYE866_Z0_DW_RATIO", "Z0")
    DYE866.write_new_commondata(
        Path("data_reimplemented_PDXSECRATIO.yaml"),
        Path("kinematics_reimplemented_PDXSECRATIO.yaml"),
        Path("uncertainties_reimplemented_PDXSECRATIO.yaml"),
    )
    if covmat_is_close("DYE866_Z0_800GEV_DW_RATIO_PDXSECRATIO", "legacy", "reimplemented"):
        print("covmat is close")
    else:
        print("covmat is different.")


if __name__ == "__main__":
    main()
