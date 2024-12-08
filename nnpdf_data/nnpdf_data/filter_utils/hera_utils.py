from pathlib import Path
from dataclasses import dataclass
import typing
from typing import List
import numpy as np
import pandas as pd
from os import PathLike
import yaml
#from validphys.api import API

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
      

# Subroutines for testing the implementation of the commondata.

## Obtain the covariance matrix for a given variant.
#def _covmat(name: str, var: str):
#   inp = {"dataset_input": {"dataset": name, "variant": var}, "theoryid": 700, "use_cuts": "internal"}
#   return API.covmat_from_systematics(**inp)
#
## Compare the covariance matrices of two different variants. True if close.
#def covmat_is_close(name: str, var1: str, var2: str) -> bool:
#   return np.isclose(_covmat(name,var1),_covmat(name,var2)).all()




