from nnpdf_data.filter_utils.hera_utils import commondata, covmat_is_close
from pathlib import Path
from dataclasses import dataclass
import typing
from typing import List
import numpy as np
import pandas as pd
from os import PathLike
from fortranformat import FortranRecordWriter
import yaml

@dataclass
class hera_commondata(commondata):
   def __init__(self, filename: str | PathLike, dataset_name: str, 
                process: str):
      # Read the data.
      file = Path(filename)
      df = pd.read_table(file, sep=r"\s+")

      # Kinematic quantieties.
      self.central_values = df["Sigma"].to_numpy()
      self.kinematics = df[["x", "Q2", "y"]].to_numpy()
      self.kinematic_quantities = ["x", "Q2", "y"]

      # Statistical uncertainties.
      statistical_uncertainties = df["stat"].to_numpy()
      for iunc,unc in enumerate(statistical_uncertainties):
         unc = self.central_values[iunc]*unc/100
         statistical_uncertainties[iunc] = unc
      self.statistical_uncertainties = statistical_uncertainties

      # Systematic uncertainties.
      # remove the column containing the total uncertainty excluding 
      # procedural uncertainties.
      df = df.drop(columns=["tot_noproc"])
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
            systypes.append(("MULT", f"HC_{name}"))
      self.systypes = systypes
      self.process = process
      self.dataset_name = dataset_name

def main():
   print("Reimplementing the HERA commondata")   
   hera_em = hera_commondata("./rawdata/HERA1+2_CCem.dat","HERACOMBCCEM", "DIS_CCE")
   hera_em.write_new_commondata(Path("data_reimplemented_EM-SIGMARED.yaml"),
                                Path("kinematics_reimplemented_EM-SIGMARED.yaml"),
                                Path("uncertainties_reimplemented_EM-SIGMARED.yaml"))
   hera_ep = hera_commondata("./rawdata/HERA1+2_CCep.dat","HERACOMBCCEP", "DIS_CCE")
   hera_ep.write_new_commondata(Path("data_reimplemented_EP-SIGMARED.yaml"),
                                Path("kinematics_reimplemented_EP-SIGMARED.yaml"),
                                Path("uncertainties_reimplemented_EP-SIGMARED.yaml"))
   # Check if the covariance matrix of the reimplemented data is close to the
   # legacy implementation
   print("Check covariance matrix for HERA_CC_318GEV_EM-SIGMARED:")
   if(covmat_is_close("HERA_CC_318GEV_EM-SIGMARED","reimplemented","legacy")):
      print("Covmat is close.")
   else:
      print("Covmat is different.")
   print("Check covariance matrix for HERA_CC_318GEV_EP-SIGMARED:")
   if(covmat_is_close("HERA_CC_318GEV_EP-SIGMARED","reimplemented","legacy")):
      print("Covmat is close.")
   else:
      print("Covmat is different.")

if __name__ == "__main__":
   main()



