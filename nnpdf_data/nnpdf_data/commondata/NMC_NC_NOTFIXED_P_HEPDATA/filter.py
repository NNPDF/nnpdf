r"""Implement data form Hepdata reference, data from Table 5.
R is simultaneously determined from the data.
Results are averaged onver different $\sqrt{s}$.
"""

import pathlib

import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent


def read_tables(store_path, header_line):
    """Parse Tables."""
    dfs = pd.DataFrame()
    for file in store_path.iterdir():
        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            df = pd.DataFrame(
                [l.split(",") for l in lines[header_line:-1]],
                columns=[
                    "Q2",
                    "R",
                    "F2",
                    "stat+",
                    "stat-",
                    "sys+",
                    "sys-",
                ],
            )
            df["x"] = float(lines[header_line - 2].split(",")[1])
            dfs = pd.concat([dfs, df], ignore_index=True) if not dfs.empty else df

    dfs = dfs.astype(float)
    return dfs.sort_values(["x", "Q2"])


def write_files(df, store_path):
    """Write kinematics, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [float(x) for x in df["F2"]]}
    with open(store_path / "data.yaml", "w", encoding="utf-8") as file:
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
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(store_path / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

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
    }
    error = []
    for _, row in df.iterrows():
        e = {
            "stat": float(row["stat+"]),
            "sys": float(row["sys+"]),
        }
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(store_path / "uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    df = read_tables(HERE / "rawdata", header_line=14)
    write_files(df, HERE)
