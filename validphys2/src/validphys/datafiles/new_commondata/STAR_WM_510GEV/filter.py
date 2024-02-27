import pandas as pd
import yaml
import glob

ECM = 510
MW = 80.398


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        df = pd.read_csv(fname, delimiter=",", skiprows=10)
        df["M2"] = MW**2
        df["sqrts"] = ECM
    return df


def write_data(df):
    data_central = df["$A_L$"].tolist()
    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w", encoding="utf-8") as file:
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

    with open("kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(len(df)):
        # here uncertainties are symmetric
        e = {
            "stat": float(df.loc[i, "stat +"]),
            "sys": float(df.loc[i, "syst +"]),
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
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    fnames = glob.glob("rawdata/*.csv")
    df = read_data(fnames)
    write_data(df)
