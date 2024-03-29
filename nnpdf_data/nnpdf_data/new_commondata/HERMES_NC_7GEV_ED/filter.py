import glob
import pathlib
import sys

import numpy as np
import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent
sys.path = [str(HERE.parent / "HERMES_NC_7GEV_EP")] + sys.path

from filter import compute_covmat


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            data = yaml.safe_load(file)

        xsub = data["independent_variables"][0]["values"]
        Qsub = data["independent_variables"][1]["values"]
        Gsub = data["dependent_variables"][1]["values"]

        for i in range(len(xsub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "x": [xsub[i]["value"]],
                            "Q2": [Qsub[i]["value"]],
                            "G": [Gsub[i]["value"]],
                            "stat": [Gsub[i]["errors"][0]["symerror"]],
                            "exp": [Gsub[i]["errors"][1]["symerror"]],
                            "param": [Gsub[i]["errors"][2]["symerror"]],
                            "evol": [Gsub[i]["errors"][3]["symerror"]],
                        }
                    ),
                ],
                ignore_index=True,
            )

    return df


def read_corrmatrix(nb_datapoints: int = 15) -> np.ndarray:
    """Load the correlation Matrix in Table 23."""
    file = pathlib.Path("./rawdata/HEPData-ins726689-v1-Table_23.yaml")
    loaded_file = yaml.safe_load(file.read_text())

    corrs = loaded_file['dependent_variables'][0]['values']
    df_corrs = pd.DataFrame(corrs)

    return df_corrs.value.values.reshape((nb_datapoints, nb_datapoints))


def write_data(df):
    data_central = []
    for i in range(len(df["G"])):
        data_central.append(float(df.loc[i, "G"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df["G"])):
        kin_value = {
            "x": {"min": None, "mid": float(df.loc[i, "x"]), "max": None},
            "Q2": {"min": None, "mid": float(df.loc[i, "Q2"]), "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Extract the correlation matrix and compute artificial systematics
    ndata_points = len(data_central)
    corrmatrix = read_corrmatrix(nb_datapoints=ndata_points)
    # Compute the covariance matrix
    compute_covmat(corrmatrix, df, ndata_points)

    # Compute the covariance matrix
    art_sys = compute_covmat(corrmatrix, df, ndata_points)

    error = []
    for i in range(ndata_points):
        e = {}
        # add the art sys
        for j in range(ndata_points):
            e[f"sys_{j}"] = art_sys[i][j]

        e["stat"] = 0  # This is set to 0 as the stat unc is correlated and reported in sys_0
        e["exp"] = float(df.loc[i, "exp"])  # experimental including normalization
        e["param"] = float(df.loc[i, "param"])
        e["evol"] = float(df.loc[i, "evol"])
        error.append(e)

    error_definition = {
        f"sys_{i}": {
            "description": f"{i} artificial correlated statistical uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(ndata_points)
    }

    error_definition.update(
        {
            "stat": {
                "description": "statistical uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "exp": {
                "description": "experimental systematic uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "param": {
                "description": "parametrization systematic uncertainty",
                "treatment": "ADD",
                "type": "CORR",
            },
            "evol": {
                "description": "evolution systematic uncertainty",
                "treatment": "ADD",
                "type": "CORR",
            },
        }
    )

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    fnames = glob.glob("./rawdata/HEPData-ins726689-v1-Table_13.yaml")
    nnames = sorted([i for i in fnames])
    df = read_data(nnames)
    write_data(df)
