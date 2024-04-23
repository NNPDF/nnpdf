import glob

import pandas as pd
import yaml
import sys

sys.path.append('../../')
from filter_utils.correlations import compute_covmat


def read_data(fnames):
    df3 = pd.DataFrame()
    df4 = pd.DataFrame()
    for fname in fnames:
        if "3" in fname:
            dfi = pd.read_csv(fname)
            if "pT" in fname:
                correlation = False
                df3["pT"] = dfi["Parton Jet $p_T$ [GeV/c]"]

            elif "ALL" in fname:
                df3["ALL"] = dfi["$A_{LL}$"]
                df3["stat_min"] = dfi["stat +"]
                df3["stat_max"] = dfi["stat -"]
                df3["sys_min"] = dfi["sys +"]
                df3["sys_max"] = dfi["sys -"]

            df3["abs_eta_min"] = 0.0
            df3["abs_eta_max"] = 0.5
            df3["sqrts"] = 200

        elif "4" in fname:
            dfi = pd.read_csv(fname)
            if "pT" in fname:
                correlation = False
                df4["pT"] = dfi["Parton Jet $p_T$ [GeV/c]"]

            if "ALL" in fname:
                df4["ALL"] = dfi["$A_{LL}$"]
                df4["stat_min"] = dfi["stat +"]
                df4["stat_max"] = dfi["stat -"]
                df4["sys_min"] = dfi["sys +"]
                df4["sys_max"] = dfi["sys -"]

            df4["abs_eta_min"] = 0.5
            df4["abs_eta_max"] = 1.0
            df4["sqrts"] = 200

        elif "5" in fname:
            correlation = True
            dfc_col = pd.read_csv(fname)

            dfc = pd.DataFrame()
            biny = 1
            for i in range(len(dfc_col)):
                if dfc_col.loc[i, "binx"] == "-":
                    biny = float(dfc_col.loc[i, "val"])
                else:
                    binx = float(dfc_col.loc[i, "binx"])
                    dfc.loc[binx, biny] = dfc_col.loc[i, "val"]

            dfc = dfc.astype("float")

        else:
            print("ERROR: Unknown table number detected! Check input files.")

    df = pd.concat([df3, df4], ignore_index=True)

    for i in range(len(df)):
        df.loc[i, "stat"] = df.loc[i, "stat_max"]
        df.loc[i, "sys"] = df.loc[i, "sys_max"]
    return df, dfc


def write_data(df, dfc):
    data_central = []
    for i in range(len(df["ALL"])):
        data_central.append(float(df.loc[i, "ALL"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Compute the covariance matrix
    corrmat = dfc.values.reshape((len(df), len(df)))
    art_sys = compute_covmat(corrmat, df["stat"], len(df))

    # Write kin file
    kin = []
    for i in range(len(df["ALL"])):
        kin_value = {
            "pT": {"min": None, "mid": float(df.loc[i, "pT"]), "max": None},
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "abs_eta": {
                "min": float(df.loc[i, "abs_eta_min"]),
                "mid": None,
                "max": float(df.loc[i, "abs_eta_max"]),
            },
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(len(df)):
        e = {}
        for j in range(len(df)):
            e[f"sys_{j}"] = art_sys[i][j]
        error.append(e)

    error_definition = {
        f"sys_{i}": {
            "description": f"{i} artificial correlated statistical uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(len(df))
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # TODO: Need to generate `observable` cards and corresponding
    # pineappl grids and FK tables as the orders have changed!!!!
    fnames = glob.glob("rawdata/*.csv")
    nnames = sorted([i for i in fnames])
    df, dfc = read_data(nnames)
    write_data(df, dfc)
