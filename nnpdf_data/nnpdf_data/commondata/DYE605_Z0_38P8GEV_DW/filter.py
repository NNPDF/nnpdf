from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import typing
from typing import List

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.hera_utils import commondata, covmat_is_close


def mergetables() -> pd.DataFrame:

    table_paths = []
    for i in range(1, 8):
        table_paths.append(Path(f"./rawdata/Table{i}.csv"))

    # List with the rapidity bins for tables 1 to 7.
    yrap = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]

    col_names = ["M2", "dsig", "statp", "statm", "normp", "normm", "sysp", "sysm"]
    col_names_all = col_names + ["y", "sqrts"]

    combined_df = pd.DataFrame(columns=col_names_all)
    for i, path in enumerate(table_paths):
        df = pd.read_csv(path, header=11, names=col_names)
        df["y"] = yrap[i]
        df["sqrts"] = 38.8
        df = df[pd.to_numeric(df['dsig'], errors='coerce').notnull()]
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    # In the table we have sqrt(tau) not M2; compute M2=tau*s
    combined_df["M2"] = (combined_df["M2"] * 38.8) ** 2

    return combined_df


def nuclear_uncert_dw(tableN: PathLike, tablep: PathLike):
    dfN = pd.read_table(tableN)
    dfp = pd.read_table(tablep)
    return dfN, dfp


@dataclass
class E605_commondata(commondata):
    def __init__(self, data: pd.DataFrame, dataset_name: str, process: str):

        # Kinematic quantities.
        self.central_values = data["dsig"].astype(float).to_numpy()
        self.kinematics = data[["y", "M2", "sqrts"]].astype(float).to_numpy()
        self.kinematic_quantities = ["y", "M2", "sqrts"]

        # Statistical uncertainties.
        self.statistical_uncertainties = data["statp"]

        # the overall 10% statistical uncertainty is treated as
        # additive, while normalisation uncertainty is always treated
        # multiplicatively
        syst = pd.DataFrame(0.1 * self.central_values)

        # Systematic uncertainties.
        syst["norm"] = self.central_values * data["normp"].str.strip("%").astype(float) / 100

        # self.systematic_uncertainties = np.dstack((stat,norm))[0]
        self.systypes = [("ADD", "UNCORR"), ("MULT", "CORR")]

        # Compute the point-to-point uncertainties
        nrep = 999
        norm = np.sqrt(nrep)
        dfN, dfp = nuclear_uncert_dw(
            "rawdata/nuclear/output/tables/group_result_table.csv",
            "rawdata/proton_ite/output/tables/group_result_table.csv",
        )

        for rep in range(1, nrep + 1):
            Delta = (dfN[f"rep_{rep:05d}"] - dfp["theory_central"]) / norm
            syst[f"NUCLEAR{rep:05d}"] = Delta
            self.systypes.append(("ADD", f"NUCLEAR{rep:05d}"))

        self.systematic_uncertainties = syst.to_numpy()

        self.process = process
        self.dataset_name = dataset_name


def main():
    data = mergetables()
    # First create the commondata variant without the nuclear uncertainties.
    DYE605 = E605_commondata(data, "DYE605_Z0_38P8GEV", "Z0")
    DYE605.write_new_commondata(
        Path("data_reimplemented_PXSEC.yaml"),
        Path("kinematics_reimplemented_PXSEC.yaml"),
        Path("uncertainties_reimplemented_PXSEC.yaml"),
    )
    if covmat_is_close("DYE605_Z0_38P8GEV_DW_PXSEC", "legacy", "reimplemented"):
        print("covmat is close")
    else:
        print("covmat is different.")


if __name__ == "__main__":
    main()
