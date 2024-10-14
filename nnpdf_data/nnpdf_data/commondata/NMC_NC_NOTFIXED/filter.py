r"""Implement data form Hepdata reference, data from Table 3."""

import pathlib

import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent


def read_tables():
    """Parse Tables."""
    dfs = pd.DataFrame()
    for table_id in range(2, 22):
        file = HERE / "rawdata" / f"Table{table_id}.csv"
        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            df = pd.DataFrame(
                [l.split(",") for l in lines[14:-1]],
                columns=[
                    "Q2",
                    "F2_ratio",
                    "stat+",
                    "stat-",
                    "sys+",
                    "sys-",
                ],
            )
            df["x"] = float(lines[12].split(",")[1])
            df["sqrts"] = float(lines[11].split(",")[1])
            dfs = pd.concat([dfs, df], ignore_index=True) if not dfs.empty else df

    dfs = dfs.astype(float)

    # dfs["y"] = dfs.Q2 /( dfs.x * dfs.sqrts**2)
    return dfs


def write_files(df):
    """Write kinematics, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [float(x) for x in df["F2_ratio"]]}
    with open(HERE / f"data.yaml", "w", encoding="utf-8") as file:
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
            "sqrts": {
                "min": None,
                "mid": float(row.sqrts),
                "max": None,
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(HERE / "kinematics.yaml", "w", encoding="utf-8") as file:
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
    with open(HERE / "uncertainties_hepdata.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    df = read_tables()
    write_files(df)
