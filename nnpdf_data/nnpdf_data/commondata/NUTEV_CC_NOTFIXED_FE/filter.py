"""
This filter is to construct an updated version of the legacy variant, starting from the old
Buildmaster rawdata which do not include nuclear uncertainties, but have the correct updated BR.
"""
import pathlib

import numpy as np
import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent

M_NEUTRON = 939.565346 * 0.001
M_PROTON = 938.272013 * 0.001
# The old implementation used the isoscalar average.
MN = (M_PROTON + M_NEUTRON) / 2


# The value of BrC 0.086 +- 0.005 according to the PDG(2017)
# The uncertainty on the BR is now accounted for as
# an additional fully correlated systematic uncertainty
BRC = 0.087
BRC_err = 0.005

# 2.1 % Normalization error, from old implementation
NORM_err = 0.021


def read_table(projectile):
    """Parse the old tables."""
    with open(HERE / "rawdata" / f"NuTeVtable_{projectile}.dat", encoding="utf-8") as f:
        raw_data = np.loadtxt(f)

    with open(
        HERE / "rawdata" / f"nf20-1.25-0.60.{projectile}.cor", encoding="utf-8"
    ) as f:
        acc_corr = np.loadtxt(f)

    df = pd.DataFrame()
    df["E"] = raw_data[:, 1]
    df["y"] = raw_data[:, 2]
    df["x"] = raw_data[:, 3]
    df["Q2"] = 2 * df.x * df.y * df.E * MN
    df["sigma"] = raw_data[:, 4]
    df["stat"] = raw_data[:, 5]
    df["sys"] = raw_data[:, 6]
    df["acc"] = acc_corr[:, 0]
    # df["acc_err"] = acc_corr[:,1] # TODO: why not used ??
    df["BRC_err"] = BRC_err / BRC * df["sigma"]
    df["norm"] = NORM_err * df["sigma"]

    return df


def write_files(df, variant):
    """Write kinematics, central value and uncertainties files."""

    NORM_FACT = 1 / (BRC * df["acc"])

    # Write central data
    data_central_yaml = {
        "data_central": [round(float(x), 11) for x in df["sigma"] * NORM_FACT]
    }
    with open(HERE / f"data_{variant}.yaml", "w", encoding="utf-8") as file:
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
                "min": None,
                "mid": float(row.y),
                "max": None,
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(HERE / f"kinematics_{variant}.yaml", "w", encoding="utf-8") as file:
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
            "description": "Normalization uncertainty",
            "treatment": "MULT",
            "type": "NUTEVNORM1",
        },
        "brunc": {
            "description": "Branching ratio uncertainty",
            "treatment": "MULT",
            "type": "NUTEVBRC1",
        },
    }
    error = []
    for k, row in df.iterrows():
        e = {
            "stat": float(row["stat"] * NORM_FACT[k]),
            "sys": float(row["sys"] * NORM_FACT[k]),
            "norm": float(row["norm"] * NORM_FACT[k]),
            "brunc": float(row["BRC_err"] * NORM_FACT[k]),
        }
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(
        HERE / f"uncertainties_{variant}_hepdata.yaml", "w", encoding="utf-8"
    ) as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    for proj, variant in [("nu", "NU"), ("bar", "NB")]:
        df = read_table(proj)
        write_files(df, variant)
