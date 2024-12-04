"""Implement Hepdata F2 charm data from Table 1 (ie. fig 14 of the paper). EMC assumes R=0.
Branching Ration of 0.82 has been included to match the normalization of the legacy version."""

import pandas as pd
import pathlib
import yaml

HERE = pathlib.Path(__file__).parent

# Normalization factor, to match old implementation
# Most likely due to the Branching Ratio
NORM = 0.82
BR_ERR = 0.012

M_P = 0.938

# NOTE: Systematic uncertainties which is estimated as in the legacy variant,
# as a factor of 15 %, see old buildmaster.
SYST_ERR = 0.15

def read_tables():
    """Read Hepdata table."""

    dfs = pd.DataFrame()
    with open(HERE / "rawdata" / "HEPData-ins180921-v1-Table_1.csv", encoding="utf-8") as f:
        lines = f.readlines()

        # loop on columns
        for idx_line in range(11, 137, 14):
            df = pd.DataFrame(
                [l.split(",") for l in lines[idx_line + 3 : idx_line + 11]],
                columns=["Q2", "Q2_min", "Q2_max", "F2", "error+", "error-"],
            )
            df = df[~df.F2.str.contains("-")].astype(float)
            df["x_min"] = float(lines[idx_line + 1].split("TO")[0].split(",")[-1])
            df["x_max"] = float(lines[idx_line + 1].split("TO")[-1])
            df["x"] = (df.x_min + df.x_max) / 2
            df["sqrts"] = float(lines[idx_line].split(",")[-1])
            df["y"] = df.Q2 / (df.x * (df.sqrts**2 - M_P**2))
            dfs = pd.concat([dfs, df], ignore_index=True) if not dfs.empty else df
    return dfs


def write_files(df):
    """Write kinematics, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [round(float(x) * NORM, 7) for x in df["F2"]]}
    with open(HERE / "data_rzero.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for _, row in df.iterrows():
        kin_value = {
            "x": {"min": float(row.x_min), "mid": float(row.x), "max": float(row.x_max)},
            "Q2": {
                "min": float(row.Q2_min),
                "mid": round(float(row.Q2), 2),
                "max": float(row.Q2_max),
            },
            "y": {"min": None, "mid": float(row.y), "max": None},
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(HERE / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # loop on data points
    error_definition = {
        "stat": {
            "description": "Total statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "sys": {
            "description": "Total systematical uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
        "norm": {
            "description": "Normalization error due to BR ratio",
            "treatment": "MULT",
            "type": "CORR",
        },
    }
    error = []
    for _, row in df.iterrows():
        e = {
            "stat": round(float(row["error+"]) * NORM, 7),
            "sys": round(float(row["F2"]) * NORM * SYST_ERR, 7),
            "norm": round(float(row["F2"]) * NORM * BR_ERR, 7),
        }
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(HERE / "uncertainties_rzero.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    df = read_tables()
    write_files(df)
