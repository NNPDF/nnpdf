"""This script provides the filer to the inclusive jet STAR 2004 dataset."""
import pandas as pd
import yaml

POL_UNC = 0.25
SQRTS = 200
ETA_MIN = 0.2
ETA_MAX = 0.8


def read_data():
    df = pd.DataFrame()
    file_raw = "HEPData-ins723509-v1-Figure_3.csv"
    with open(f"rawdata/{file_raw}", "r", encoding="utf-8") as file:
        df = pd.read_csv(file, skiprows=8)

    df["stat"] = df["stat +"]
    df["sys"] = df["syst +"]
    df["pol"] = POL_UNC * abs(df["$A_{LL}$"])
    df["sqrts"] = SQRTS
    df["eta"] = (ETA_MAX + ETA_MIN) / 2
    return df


def write_data(df):
    data_central = []
    for i in range(len(df["$A_{LL}$"])):
        data_central.append(float(df.loc[i, "$A_{LL}$"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df["$A_{LL}$"])):
        kin_value = {
            "pT": {
                "min": float(df.loc[i, "$p_T (GeV/c)$ LOW"]),
                "mid": float(df.loc[i, "$p_T (GeV/c)$"]),
                "max": float(df.loc[i, "$p_T (GeV/c)$ HIGH"]),
            },
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "eta": {
                "min": ETA_MIN,
                "mid": float(df.loc[i, "eta"]),
                "max": ETA_MAX,
            },
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w", encoding="utf-8") as file:
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
            "type": "RHIC2004POL",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    df = read_data()
    write_data(df)
