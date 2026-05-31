from dataclasses import dataclass
from os import PathLike
from pathlib import Path

import pandas as pd

from nnpdf_data.filter_utils.hera_utils import commondata  # , covmat_is_close


@dataclass
class hera_commondata(commondata):
    def __init__(self, filename: str | PathLike, dataset_name: str, process: str):
        # Read the data.
        file = Path(filename)
        df = pd.read_table(file, sep=r"\s+")

        # Kinematic quantieties.
        self.central_values = df["Sigma"].to_numpy()
        self.kinematics = df[["x", "Q2", "y"]].to_numpy()
        self.kinematic_quantities = ["x", "Q2", "y"]

        # Statistical uncertainties.
        statistical_uncertainties = df["stat"].to_numpy(copy=True)
        for iunc, unc in enumerate(statistical_uncertainties):
            unc = self.central_values[iunc] * unc / 100
            statistical_uncertainties[iunc] = unc
        self.statistical_uncertainties = statistical_uncertainties

        # Systematic uncertainties.
        # remove the column containing the total uncertainty excluding
        # procedural uncertainties.
        df = df.drop(columns=["tot_noproc"])
        sys_uncert_col_names = list(df.columns.values)[5:]
        self.systematic_uncertainties = df[sys_uncert_col_names].to_numpy()
        systematic_uncertainties = df[sys_uncert_col_names].to_numpy()
        for iunc, unc in enumerate(systematic_uncertainties):
            unc = self.central_values[iunc] * unc / 100
            systematic_uncertainties[iunc] = unc
        self.systematic_uncertainties = systematic_uncertainties

        # All uncertainties are treated as multiplicative.
        systypes = []
        for name in sys_uncert_col_names:
            if name == "uncor":
                systypes.append(("MULT", "UNCORR"))
            else:
                systypes.append(("MULT", "HC_" + name))
        self.systypes = systypes
        self.process = process
        self.dataset_name = dataset_name


def main():
    hera_ep = hera_commondata("./rawdata/HERA1+2_NCep_460.dat", "HERACOMBNCEP460", "DIS_NCE")
    hera_ep.write_new_commondata(
        Path("data_EP-SIGMARED.yaml"),
        Path("kinematics_EP-SIGMARED.yaml"),
        Path("uncertainties_EP-SIGMARED.yaml"),
    )


if __name__ == "__main__":
    main()
