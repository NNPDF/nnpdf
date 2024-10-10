"""Implement BCDMS_NC_NOTFIXED_P data from hepdata. We use data with R(QCD)."""

import pathlib

import numpy as np
import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent
VARIANT = "hepdata"


def read_tables():
    """Parse Tables."""
    dfs = pd.DataFrame()
    for file in pathlib.Path(HERE / "rawdata").iterdir():
        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            df = pd.DataFrame(
                [l.split(",") for l in lines[14:-1]],
                columns=[
                    "Q2",
                    "F2",
                    "stat+",
                    "stat-",
                    "sys+",
                    "sys-",
                    "norm+",
                    "norm-",
                ],
            )
            df["x"] = float(lines[12].split(",")[1])
            try:
                df["sqrts"] = float(lines[11].split(",")[1])
                df["sqrts_min"] = df["sqrts"]
                df["sqrts_max"] = df["sqrts"]
            except ValueError:
                df["sqrts_min"] = float(lines[11].split(",")[1].split("-")[0])
                df["sqrts_max"] = float(lines[11].split(",")[1].split("-")[1])
                df["sqrts"] = (df["sqrts_min"] + df["sqrts_max"]) / 2
            if dfs.empty:
                dfs = df
            else:
                dfs = pd.concat([dfs, df], ignore_index=True)
    dfs["norm+"] = [float(x[:-1]) for x in dfs["norm+"]]
    dfs["norm-"] = [float(x[:-2]) for x in dfs["norm-"]]
    dfs = dfs.astype(float)
    dfs["norm+"] *= abs(dfs.F2)
    dfs["norm-"] *= abs(dfs.F2)

    # dfs["y"] = dfs.Q2 /( dfs.x * dfs.sqrts**2) 
    return dfs.sort_values(["Q2", "x", "sqrts"])


def write_files(df):
    """Write kinemati, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [float(x) for x in df["F2"]]}
    with open(HERE / f"data_{VARIANT}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for _, row in df.iterrows():
        kin_value = {
            "x": {
                "min": None,
                "mid": float(row.x),
                "max": None,
            },
            "Q2": {
                "min": None,
                "mid": float(row.Q2),
                "max": None,
            },
            "y": {
                "min": float(row.sqrts_min),
                "mid": float(row.sqrts),
                "max": float(row.sqrts_max),
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(HERE / f"kinematics_{VARIANT}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file)

    # loop on data points
    error_definition = {
        "stat": {
            "description": "statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "sys": {
            "description": "systematical uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
        "norm": {
            "description": "normalization uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
    }
    error = []
    for _, row in df.iterrows():
        e = {
            "stat": float(row["stat+"]),
            "sys": float(row["sys+"]),
            "norm": float(row["norm+"]),
        }
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(HERE / f"uncertainties_{VARIANT}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    df = read_tables()
    write_files(df)
