from pathlib import Path

import pandas as pd
import yaml


def read_data(path_rawdata: str) -> pd.DataFrame:
    data = yaml.safe_load(Path(path_rawdata).read_text())

    ptvals = data["independent_variables"][0]["values"]
    asym_obs = data["dependent_variables"][0]["values"]

    concatenated_table = []
    for i in range(len(ptvals)):
        # Compute pT mid-value if not given
        pt_mid = (float(ptvals[i]["high"]) + float(ptvals[i]["low"])) / 2
        pt_midval = ptvals[i].get("value", pt_mid)

        concatenated_table.append(
            pd.DataFrame(
                {
                    "pT_low": [ptvals[i]["low"]],
                    "pT": [pt_midval],
                    "pT_high": [ptvals[i]["high"]],
                    "asym": [asym_obs[i]["value"]],
                    "stat_err": [asym_obs[i]["errors"][0]["symerror"]],
                    "lumi_err": [asym_obs[i]["errors"][1]["symerror"]],
                    "spol_err": [asym_obs[i]["errors"][2]["symerror"].removesuffix("%")],
                }
            )
        )

    return pd.concat(concatenated_table, ignore_index=True)


def dump_data(df_table: pd.DataFrame) -> None:
    # Dump central data into Yaml file
    data_central = []
    for i in range(len(df_table["asym"])):
        data_central.append(float(df_table.loc[i, "asym"]))

    with open("data.yaml", "w") as file:
        yaml.dump({"data_central": data_central}, file, sort_keys=False)

    # Dump the kinematics into Yaml file
    kinematics = []
    for i in range(len(df_table["asym"])):
        kin_value = {
            "pT": {
                "min": float(df_table.loc[i, "pT_low"]),
                "mid": float(df_table.loc[i, "pT"]),
                "max": float(df_table.loc[i, "pT_high"]),
            },
            "sqrts": {"min": None, "mid": 200.0, "max": None},
        }
        kinematics.append(kin_value)

    with open("kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    # Dump the uncertainties into Yaml file
    errors = []
    for i in range(len(df_table)):
        error_per_bin = {
            "stat": float(df_table.loc[i, "stat_err"]),
            "sys_lumi": float(df_table.loc[i, "lumi_err"]),
            "sys_pol": data_central[i] * float(df_table.loc[i, "spol_err"]) / 100.0,
        }
        errors.append(error_per_bin)

    error_definition = {
        "stat": {"description": "Statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "sys_lumi": {
            "description": "Systematic uncertainties due to luminosity",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_pol": {
            "description": "Systematic uncertainties due to polarization",
            "treatment": "MULT",
            "type": "CORR",
        },
    }

    with open("uncertainties.yaml", "w") as file:
        yaml.dump({"definitions": error_definition, "bins": errors}, file, sort_keys=False)

    return


if __name__ == "__main__":
    df_table = read_data("./rawdata/HEPData-ins1282448-v1-Table_7.yaml")
    dump_data(df_table=df_table)
