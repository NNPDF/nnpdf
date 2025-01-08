from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import typing
from typing import List

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.hera_utils import commondata, covmat_is_close

# TODO the same function might be placed differently depending on the version of the code.
try:
    from validphys.newcommondata_utils import covmat_to_artunc
except ModuleNotFoundError:
    from nnpdf_data.filter_utils.correlations import covmat_to_artunc


def readdata() -> pd.DataFrame:
    col_names = ["xtlow", "xtup", "xt", "xb", "M", "pt", "sigr", "stat", "sys", "dx"]
    table_path = Path(f"./rawdata/table.csv")
    df = pd.read_csv(table_path, names=col_names, header=0)
    return df


def nuclear_uncert_dw(tableN: PathLike, tablep: PathLike):
    dfN = pd.read_table(tableN)
    dfp = pd.read_table(tablep)
    return dfN, dfp


def get_covmat_list(covmat_filename: PathLike):
    covmat = np.loadtxt(covmat_filename, delimiter=",")
    return covmat


def is_symmetric(mat: np.ndarray) -> bool:
    return np.all(np.isclose(mat, mat.T))


@dataclass
class E906_DW_RATIO_commondata(commondata):
    def __init__(self, data: pd.DataFrame, dataset_name: str, process: str):

        self.central_values = data["sigr"].astype(float).to_numpy()
        # Kinematic quantities.
        y = 0.5 * np.log(data["xb"] / data["xt"])
        Ebeam = 120
        MP = 0.938
        s = 2 * MP * (Ebeam + MP) * y / y

        kin = pd.concat([y, data["M"] ** 2, s], axis=1)
        kin = kin.set_axis(["y", "M2", "s"], axis=1)
        self.kinematics = kin.astype(float).to_numpy()
        self.kinematic_quantities = ["y", "M2", "sqrts"]

        # Statistical uncertainties.
        self.statistical_uncertainties = np.zeros(len(data["stat"]))

        # Systematic uncertainties.
        syst = pd.DataFrame(data["sys"])
        syst = syst.set_axis(["sys_unc_0"], axis=1)
        systypes = [("ADD", "CORR")]

        # Compute the artificial uncertainties from the
        # covariance matrix found in equation 9 of
        # https://arxiv.org/abs/2103.04024
        covmat = get_covmat_list("rawdata/covmat.csv")
        if not is_symmetric(covmat):
            raise ValueError("The covariance matrix is not symmetric.")
        artunc = pd.DataFrame(covmat_to_artunc(len(covmat), covmat.flatten().tolist()))
        artunc_names = []
        for i in range(1, len(artunc) + 1):
            artunc_names.append(f"sys_unc_{i}")
        artunc = artunc.set_axis(artunc_names, axis=1)

        for _, col_name in enumerate(artunc):
            systypes.append(("ADD", "CORR"))
            syst[col_name] = artunc[col_name]

        # Compute the point-to-point uncertainties
        nrep = 100
        norm = np.sqrt(nrep)
        norm = 100
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
    DYE906 = E906_DW_RATIO_commondata(data, "DYE906_Z0", "Z0")
    DYE906.write_new_commondata(
        Path("data_reimplemented_PDXSECRATIO.yaml"),
        Path("kinematics_reimplemented_PDXSECRATIO.yaml"),
        Path("uncertainties_reimplemented_PDXSECRATIO.yaml"),
    )
    if covmat_is_close("DYE906_Z0_120GEV_DW_PDXSECRATIO", "legacy", "reimplemented"):
        print("covmat is close")
    else:
        print("covmat is different.")


if __name__ == "__main__":
    main()
