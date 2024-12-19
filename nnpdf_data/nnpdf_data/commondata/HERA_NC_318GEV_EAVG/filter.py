from nnpdf_data.filter_utils.hera_utils import commondata #, covmat_is_close
from pathlib import Path
from dataclasses import dataclass
import typing
from typing import List
import numpy as np
import pandas as pd
from os import PathLike
import yaml

@dataclass
class hera_cb_commondata(commondata):
   def __init__(self, filename: str | PathLike, dataset_name: str, 
                process: str):
      # Read the data.
      file = Path(filename)
      df = pd.read_table(file, sep=r"\s+",skiprows=36)

      # Kinematic quantieties.
      self.central_values = df["Sigma"].to_numpy()
      # compute y = Q2/x/S
      S=318**2 # GeV**2
      y=df["Q2"]/df["x"]/S
      df.insert(1,"y",y.to_list())
      self.kinematics = df[["x", "Q2", "y"]].to_numpy()
      self.kinematic_quantities = ["x", "Q2", "y"]

      # Statistical uncertainties.
      statistical_uncertainties = df["stat"].to_numpy()
      for iunc,unc in enumerate(statistical_uncertainties):
         unc = self.central_values[iunc]*unc/100
         statistical_uncertainties[iunc] = unc
      self.statistical_uncertainties = statistical_uncertainties

      # Systematic uncertainties.
      sys_uncert_col_names = list(df.columns.values)[5:]
      self.systematic_uncertainties = df[sys_uncert_col_names].to_numpy()
      systematic_uncertainties = df[sys_uncert_col_names].to_numpy()
      for iunc,unc in enumerate(systematic_uncertainties):
         unc = self.central_values[iunc]*unc/100
         systematic_uncertainties[iunc] = unc
      self.systematic_uncertainties = systematic_uncertainties
      
      # All uncertainties are treated as multiplicative.
      systypes = []
      for name in sys_uncert_col_names:
         if(name == "uncor"):
            systypes.append(("MULT", "UNCORR"))
         else:
            systypes.append(("MULT", "HC_" + name))
      self.systypes = systypes
      self.process = process
      self.dataset_name = dataset_name


def main():
   hera_b = hera_cb_commondata("./rawdata/d18-037.tableBeauty.txt","HERACOMBNCEP", "DIS_NCE")
   hera_b.write_new_commondata(Path("data_BOTTOM-SIGMARED.yaml"),
                                Path("kinematics_BOTTOM-SIGMARED.yaml"),
                                Path("uncertainties_BOTTOM-SIGMARED.yaml"))

   hera_c = hera_cb_commondata("./rawdata/d18-037.tableCharm.txt","HERACOMBNCEP", "DIS_NCE")
   hera_c.write_new_commondata(Path("data_CHARM-SIGMARED.yaml"),
                                Path("kinematics_CHARM-SIGMARED.yaml"),
                                Path("uncertainties_CHARM-SIGMARED.yaml"))

if __name__ == "__main__":
   main()

