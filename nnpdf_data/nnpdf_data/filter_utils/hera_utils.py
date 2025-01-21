from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
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
    def write_new_commondata(
        self, data_filename: Path, kinematics_filename: Path, uncertainties_filename: Path
    ):
        # central data values
        data = {"data_central": self.central_values.tolist()}
        with data_filename.open("w+") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

        # kinematic quantieties
        # TODO add arrays for min and max values to derived type?
        bins = []
        for kin in self.kinematics.tolist():
            bins.append(
                {
                    self.kinematic_quantities[0]: {"min": None, "mid": kin[0], "max": None},
                    self.kinematic_quantities[1]: {"min": None, "mid": kin[1], "max": None},
                    self.kinematic_quantities[2]: {"min": None, "mid": kin[2], "max": None},
                }
            )
        data = {"bins": bins}
        with kinematics_filename.open("w+") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

        # uncertainties
        # There is only one statistical uncertainty per datapoint.
        definitions = {
            "stat": {
                "description": "Statistical uncertainty.",
                "treatment": "ADD",
                "type": "UNCORR",
            }
        }
        for isys, sys in enumerate(self.systypes):
            definitions.update(
                {
                    f"sys_corr_{isys}": {
                        "description": f"Systematic uncertainty {isys}",
                        "treatment": sys[0],
                        "type": sys[1],
                    }
                }
            )
        bins = {"bins": []}
        for i, _ in enumerate(self.central_values):
            systematics = {"stat": self.statistical_uncertainties.tolist()[i]}
            for isys, sys in enumerate(self.systematic_uncertainties[i].tolist()):
                systematics.update({f"sys_corr_{isys}": sys})
            bins["bins"].append(systematics)
        data = {"definitions": definitions}

        with uncertainties_filename.open("w+") as f:
            yaml.dump(data, f, sort_keys=False)
            yaml.dump(bins, f, sort_keys=False)
