import glob

import pandas as pd
import yaml


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            data = yaml.safe_load(file)

        pTbsub = data["independent_variables"][0]["values"]
        pTsub = data["dependent_variables"][0]["values"]
        ALLsub = data["dependent_variables"][1]["values"]

        for i in range(len(ALLsub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "pT": [pTsub[i]["value"]],
                            "pTmin": [pTbsub[i]["low"]],
                            "pTmax": [pTbsub[i]["high"]],
                            "eta": [0.0],
                            "sqrts": [200],
                            "ALL": [ALLsub[i]["value"]],
                            "stat": [ALLsub[i]["errors"][0]["symerror"]],
                        }
                    ),
                ],
                ignore_index=True,
            )

    return df


def write_data(df):
    data_central = []
    for i in range(len(df["ALL"])):
        data_central.append(float(df.loc[i, "ALL"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df["ALL"])):
        kin_value = {
            "pT": {"min": float(df.loc[i, "pTmin"]), "mid": float(df.loc[i, "pT"]), "max": float(df.loc[i, "pTmax"])},
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "eta": {"min": None, "mid": float(df.loc[i, "eta"]), "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for idx, i in enumerate(range(len(df))):
        e = {
            "stat": float(df.loc[i, "stat"]),
        }
        error.append(e)

    error_definition = {
        "stat": {"description": "statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # TODO: Need to generate `observable` cards and corresponding
    # pineappl grids and FK tables as the orders have changed!!!!
    fnames = glob.glob("rawdata/*.yaml")
    nnames = sorted([i for i in fnames])
    df = read_data(nnames)
    write_data(df)
