import glob
import sys

import pandas as pd
import yaml

from nnpdf_data.filter_utils.uncertainties import symmetrize_errors as se
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

POL_UNC = 0.094


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            data = yaml.safe_load(file)
            if "14" in fname:
                eta_min = 0.2
                eta_max = 0.8

            elif "15" in fname:
                eta_min = -0.7
                eta_max = 0.9

            else:
                print("ERROR: Unknown table number detected! Check input files.")

        pTsub = data["independent_variables"][0]["values"]
        ALLsub = data["dependent_variables"][0]["values"]

        for i in range(len(ALLsub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "pT": [pTsub[i]["value"]],
                            "pTmin": [pTsub[i]["low"]],
                            "pTmax": [pTsub[i]["high"]],
                            "eta": [(eta_min + eta_max) / 2],
                            "eta_min": [eta_min],
                            "eta_max": [eta_max],
                            "sqrts": [200],
                            "ALL": [ALLsub[i]["value"] * 1e-3],
                            "stat": [ALLsub[i]["errors"]["symerror"]["value"] * 1e-3],
                            "sys_min": [ALLsub[i]["errors"]["asymerror"]["minus"] * 1e-3],
                            "sys_max": [ALLsub[i]["errors"]["asymerror"]["plus"] * 1e-3],
                        }
                    ),
                ],
                ignore_index=True,
            )
        for i in range(len(df)):
            shift, unc = se(df.loc[i, "sys_max"], df.loc[i, "sys_min"])
            df.loc[i, "sys"] = unc
            df.loc[i, "ALL"] += shift

    df["pol"] = POL_UNC * abs(df["ALL"])
    return df


def write_data(df):
    data_central = []
    for i in range(len(df["ALL"])):
        data_central.append(float(df.loc[i, "ALL"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df["ALL"])):
        kin_value = {
            "pT": {
                "min": float(df.loc[i, "pTmin"]),
                "mid": float(df.loc[i, "pT"]),
                "max": float(df.loc[i, "pTmax"]),
            },
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "eta": {
                "min": float(df.loc[i, "eta_min"]),
                "mid": float(df.loc[i, "eta"]),
                "max": float(df.loc[i, "eta_max"]),
            },
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(len(df)):
        e = {
            "stat": float(df.loc[i, "stat"]),
            "sys": float(df.loc[i, "sys"]),
            "pol": float(df.loc[i, "pol"]),
        }
        error.append(e)

    error_definition = {
        "stat": {"description": "statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "sys": {"description": "systematic uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": "RHIC2005POL",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # TODO: Need to generate `observable` cards and corresponding
    # pineappl grids and FK tables as the orders have changed!!!!
    fnames = glob.glob("rawdata/*.yaml")
    nnames = sorted([i for i in fnames])
    df = read_data(nnames)
    write_data(df)
