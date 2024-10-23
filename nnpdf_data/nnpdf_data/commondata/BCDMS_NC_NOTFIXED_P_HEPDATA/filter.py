r"""
Implement BCDMS_NC_NOTFIXED_P_HEPDATA data form Hepdata reference. 
We use tables with R=R(QCD) and R=0, for the averaged values on $\sqrt{s}$.
The legacy implementation of BCDMS_NC_NOTFIXED_P is given by the not averaged $\sqrt{s}$ for R=0, 
so it has almost twice number of datapoints.
"""

import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import check_xq2_degenearcy

HERE = pathlib.Path(__file__).parent
VARIANTS = {"rqcd": (13, 23), "rzero": (2, 12)}


def read_tables(tables):
    """Parse Tables."""
    dfs = pd.DataFrame()
    for table_id in range(tables[0], tables[1] + 1):
        file = HERE / "rawdata" / f"Table{table_id}.csv"
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

    check_xq2_degenearcy(dfs.Q2.values, dfs.x.values)
    return dfs.sort_values(["Q2", "x", "sqrts"])


def write_files(df, variant):
    """Write kinematics, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [float(x) for x in df["F2"]]}
    with open(HERE / f"data_{variant}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for _, row in df.iterrows():
        kin_value = {
            "Q2": {
                "min": None,
                "mid": float(row.Q2),
                "max": None,
            },
            "x": {
                "min": None,
                "mid": float(row.x),
                "max": None,
            },
            "sqrts": {
                "min": float(row.sqrts_min),
                "mid": float(row.sqrts),
                "max": float(row.sqrts_max),
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
    with open(HERE / f"uncertainties_{variant}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    dfs = []
    for tables in VARIANTS.values():
        dfs.append(read_tables(tables))

    # check kinematic is the same for the 2 variants
    np.testing.assert_allclose(dfs[0].x, dfs[1].x)
    np.testing.assert_allclose(dfs[0].Q2, dfs[1].Q2)
    np.testing.assert_allclose(dfs[0].sqrts, dfs[1].sqrts)

    for df, variant in zip(dfs, VARIANTS):
        write_files(df, variant)
