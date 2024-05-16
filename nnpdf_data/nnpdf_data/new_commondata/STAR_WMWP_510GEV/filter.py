"""This script provides the common filer to the DY, WP, WM combined STAR 2011-2013 datasets.

NOTE: 
    * The beam polarization uncertainty are not included in the systematics, 
    so it is added manually.
    * Correlation are provided only for the 20213 data, so are not included.
"""
import pandas as pd
import yaml

ECM = 510
MW = 80.398
POL_UNC = 0.033


def read_data(fname):
    df = pd.read_csv(fname, delimiter=",", skiprows=10)
    df["M2"] = MW**2
    df["sqrts"] = ECM
    df["pol"] = abs(POL_UNC * df["$A_L$"])
    return df


def write_data(df, boson):
    data_central = df["$A_L$"].tolist()
    data_central_yaml = {"data_central": data_central}
    with open(f"data_{boson}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df)):
        kin_value = {
            "eta": {"min": None, "mid": float(df.loc[i, "$\eta_e$"]), "max": None},
            "M2": {"min": None, "mid": float(df.loc[i, "M2"]), "max": None},
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}

    with open(f"kinematics_{boson}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(len(df)):
        # here uncertainties are symmetric
        e = {
            "stat": float(df.loc[i, "stat +"]),
            "sys": float(df.loc[i, "syst +"]),
            "pol": float(df.loc[i, "pol"]),
        }
        error.append(e)

    error_definition = {
        "stat": {
            "description": "statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "sys": {
            "description": "systematic uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": "STARWMWPPOL",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open(f"uncertainties_{boson}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # Wp
    df = read_data("rawdata/Figure5,A_LforW^+rightarrowe^+,combineddatasamples.csv")
    write_data(df, boson="wp")
    # Wm
    df = read_data("rawdata/Figure5,A_LforW^-rightarrowe^-,combineddatasamples.csv")
    write_data(df, boson="wm")
