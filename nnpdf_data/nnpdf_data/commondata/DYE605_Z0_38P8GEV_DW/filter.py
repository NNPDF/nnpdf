from nnpdf_data.filter_utils.hera_utils import commondata #, covmat_is_close
from pathlib import Path
from dataclasses import dataclass
import typing
from typing import List
import numpy as np
import pandas as pd
from os import PathLike
import yaml

def mergetables() -> pd.DataFrame:

   table_paths = []
   for i in range(1,8):
      table_paths.append(Path(f"./rawdata/Table{i}.csv"))

   # List with the rapidity bins for tables 1 to 7.
   yrap = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]  

   col_names = ["M2","dsig","statp","statm","normp","normm","sysp","sysm"]
   col_names_all = col_names + ["y", "sqrts"]

   combined_df = pd.DataFrame(columns=col_names_all)
   for i, path in enumerate(table_paths):
      df = pd.read_csv(path, header=11, names=col_names)
      df["y"]=yrap[i]
      df["sqrts"]=38.8
      df = df[pd.to_numeric(df['dsig'], errors='coerce').notnull()]
      combined_df = pd.concat([combined_df,df],ignore_index=True)

   # In the table we have sqrt(tau) not M2; compute M2=tau*s
   combined_df["M2"] = (combined_df["M2"]*38.8)**2

   return combined_df

@dataclass
class E605_commondata(commondata):
   def __init__(self, data: pd.DataFrame, dataset_name: str, process: str):

      # Kinematic quantities.
      self.central_values = data["dsig"].astype(float).to_numpy()
      self.kinematics = data[["y", "M2", "sqrts"]].astype(float).to_numpy()
      self.kinematic_quantities = ["y", "M2", "sqrts"]

      # Statistical uncertainties.
      self.statistical_uncertainties = data["statp"]

      # Systematic uncertainties.
      norm = data["normp"].str.strip("%").astype(float).to_numpy()/100
      stat = norm/norm*0.1 # overall 10% uncertainty

      # the overall 10% statistical uncertainty is treated as
      # additive, while normalisation uncertainty is always treated
      # multiplicatively
      stat = stat * self.central_values

      self.systematic_uncertainties = np.dstack((stat,norm))[0]
      self.systypes = [("ADD", "UNCORR"),("MULT", "CORR")]

      self.process = process
      self.dataset_name = dataset_name

def main():
   data = mergetables()
   # First create the commondata variant without the nuclear uncertainties.
   DYE605 = E605_commondata(data, "DYE605_Z0_38P8GEV", "Z0")
   DYE605.write_new_commondata(Path("data_reimplemented_PXSEC.yaml"),
                                Path("kinematics_reimplemented_PXSEC.yaml"),
                                Path("uncertainties_reimplemented_PXSEC.yaml"))

if __name__ == "__main__":
   main()




