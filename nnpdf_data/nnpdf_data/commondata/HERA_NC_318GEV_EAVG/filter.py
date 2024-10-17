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
class commondata:
   central_values: np.ndarray
   kinematics: np.ndarray
   statistical_uncertainties: np.ndarray
   systematic_uncertainties: np.ndarray
   systypes: List[tuple[str, str]]
   process: str
   dataset_name: str
   kinematic_quantities: List[str]

  
   # Procedure to create data_*.yaml, kinematics_*.yaml and uncertainties_*.yaml
   def write_new_commondata(self, data_filename: str | PathLike,
                            kinematics_filename: str | PathLike,
                            uncertainties_filename: str | PathLike): 
      # central data values
      data = {"data_central": self.central_values.tolist()}
      with data_filename.open("w+") as f:
         yaml.dump(data, f, default_flow_style=False, sort_keys=False)

      # kinematic quantieties
      # TODO add arrays for min and max values to derived type?
      bins = []
      for kin in self.kinematics.tolist():
         bins.append(
            {self.kinematic_quantities[0]: 
               {
               "min": None,
               "mid": kin[0],
               "max": None
               },
            self.kinematic_quantities[1]:
               {
               "min": None,
               "mid": kin[1],
               "max": None
               },
            self.kinematic_quantities[2]:
               {
               "min": None,
               "mid": kin[2],
               "max": None
               }
             })
      data = {"bins": bins}
      with kinematics_filename.open("w+") as f:
         yaml.dump(data, f, default_flow_style=False, sort_keys=False)

      # uncertainties
      # There is only one statistical uncertainty per datapoint.
      definitions = {"stat": 
                       {
                          "description": "Statistical uncertainty.", 
                          "treatment": "ADD", 
                          "type": "UNCORR"
                       }
                    }
      for isys, sys in enumerate(self.systypes):
         definitions.update(
            {f"sys_corr_{isys}":
               {
                  "description": f"Systematic uncertainty {isys}",
                  "treatment": sys[0],
                  "type": sys[1]
               }
            })
      bins = {"bins": [] }
      for i, _ in enumerate(self.central_values):
         systematics = {"stat": self.statistical_uncertainties.tolist()[i]}
         for isys, sys in enumerate(self.systematic_uncertainties[i].tolist()):
            systematics.update({f"sys_corr_{isys}": sys})
         bins["bins"].append(systematics)
      data = {"definitions": definitions }
      # TODO Notation of reals is inconsistent from yaml.safe_dump
      #      sometimes it is in scientific notation sometimes not...
      with uncertainties_filename.open("w+") as f:
         yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)
         yaml.safe_dump(bins, f, default_flow_style=False, sort_keys=False)
      




   # TODO: old commondata format stores the uncertainties as 
   #       both additive and multiplicative. 
   def write_old_commondata(self, data_filename: str | PathLike, 
                            systype_filename: str | PathLike):
      with data_filename.open("w+") as f:
         f.write(
          f"{self.dataset_name} {len(self.systypes)} {len(self.central_values)} \n")
         FMT = "(I4,A10,"+str(3+1+1+len(self.systypes))+"E23.15)"
         print(FMT)
         line = FortranRecordWriter(FMT)
         for i in range(len(self.central_values)):
            l = ([i+1]+self.kinematics[i].tolist()+
                [self.central_values[i].tolist()]+
                [self.statistical_uncertainties[i].tolist()]+
                self.systematic_uncertainties[i].tolist())
            f.write(line.write(l)+"\n")

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
   print(" Reimplementing the HERA commondata")   
   hera_b = hera_cb_commondata("./rawdata/d18-037.tableBeauty.txt","HERACOMBNCEP", "DIS_NCE")
   hera_b.write_new_commondata(Path("data_reimplemented_BOTTOM-SIGMARED.yaml"),
                                Path("kinematics_reimplemented_BOTTOM-SIGMARED.yaml"),
                                Path("uncertainties_reimplemented_BOTTOM-SIGMARED.yaml"))

   hera_c = hera_cb_commondata("./rawdata/d18-037.tableCharm.txt","HERACOMBNCEP", "DIS_NCE")
   hera_c.write_new_commondata(Path("data_reimplemented_CHARM-SIGMARED.yaml"),
                                Path("kinematics_reimplemented_CHARM-SIGMARED.yaml"),
                                Path("uncertainties_reimplemented_CHARM-SIGMARED.yaml"))

if __name__ == "__main__":
   main()

