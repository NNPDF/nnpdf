from pathlib import Path

import pandas as pd
import yaml


def read_data(path_rawdata: str) -> pd.DataFrame:
    data = yaml.safe_load(Path(path_rawdata).read_text())

    ptvals = data["independent_variables"][0]["values"]
    asym_obs = data["dependent_variables"][0]["values"]

    concatenated_table = []
    for i in range(len(ptvals)):
        concatenated_table.append(
            pd.DataFrame(
                {
                    "pT_low": [ptvals[i]["low"]],
                    "pT_high": [ptvals[i]["high"]],
                    "asym": [asym_obs[i]["value"]],
                    "stat_err": [asym_obs[i]["errors"][0]["symerror"]],
                    "syst_err": [asym_obs[i]["errors"][1]["symerror"]],
                }
            )
        )

    return pd.concat(concatenated_table, ignore_index=True)


def dump_data(df_table: pd.DataFrame, tableid: int) -> None:
    suffix = "lowrap.yaml" if tableid == 1 else "highrap.yaml"
    # Dump central data into Yaml file
    data_central = []
    for i in range(len(df_table["asym"])):
        data_central.append(float(df_table.loc[i, "asym"]))

    with open(f"data_{suffix}", "w") as file:
        yaml.dump({"data_central": data_central}, file, sort_keys=False)

    # Dump the kinematics into Yaml file
    kinematics = []
    for i in range(len(df_table["asym"])):
        kin_value = {
            "pT": {
                "min": float(df_table.loc[i, "pT_low"]),
                "mid": (float(df_table.loc[i, "pT_high"]) + float(df_table.loc[i, "pT_low"])) / 2,
                "max": float(df_table.loc[i, "pT_high"]),
            },
            "eta": {"min": 3.15, "mid": 3.525, "max": 3.90},
        }
        kinematics.append(kin_value)

    with open(f"kinematics_{suffix}", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    # Dump the uncertainties into Yaml file
    errors = []
    for i in range(len(df_table)):
        error_per_bin = {
            "stat": float(df_table.loc[i, "stat_err"]),
            "syst": float(df_table.loc[i, "syst_err"]),
            "sys_pol": abs(data_central[i]) * 6.7 / 100.0,
        }
        errors.append(error_per_bin)

    error_definition = {
        "stat": {"description": "Statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "syst": {
            "description": "Total systematic uncertainties",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_pol": {
            "description": "Systematic uncertainties due to beam polarization",
            "treatment": "MULT",
            "type": "STAR2013POL",
        },
    }

    with open(f"uncertainties_{suffix}", "w") as file:
        yaml.dump({"definitions": error_definition, "bins": errors}, file, sort_keys=False)

    return


if __name__ == "__main__":
    for tabid in [1, 2]:
        df_table = read_data(f"./rawdata/HEPData-ins1674826-v1-Table_{tabid}.yaml")
        dump_data(df_table=df_table, tableid=tabid)
