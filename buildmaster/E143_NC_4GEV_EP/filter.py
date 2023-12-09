import pandas as pd
import yaml
import glob


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            data = yaml.safe_load(file)

        xsub = data["independent_variables"][0]["values"]
        y = 0.0
        Qsub = data["independent_variables"][1]["values"]
        Gsub = data["dependent_variables"][0]["values"]

        for i in range(len(xsub)):
            try:
                xsub[i]["low"]
            except NameError:
                xsub[i]["low"] = None
            try:
                xsub[i]["high"]
            except NameError:
                xsub[i]["high"] = None
            try:
                xsub[i]["value"]
            except KeyError:
                xsub[i]["value"] = str(
                    (float(xsub[i]["high"]) + float(xsub[i]["low"])) / 2
                )
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "x": [xsub[i]["value"]],
                            "x_low": [xsub[i]["low"]],
                            "x_high": [xsub[i]["high"]],
                            "y": [y],
                            "Q2": [Qsub[i]["value"]],
                            "G": [Gsub[i]["value"]],
                            "stat": [Gsub[i]["errors"][0]["symerror"]],
                            "sys": [Gsub[i]["errors"][1]["symerror"]],
                        }
                    ),
                ],
                ignore_index=True,
            )

    return df


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
            "x": {
                "min": float(df.loc[i, "x_low"]),
                "mid": float(df.loc[i, "x"]),
                "max": float(df.loc[i, "x_high"]),
            },
            "Q2": {"min": None, "mid": float(df.loc[i, "Q2"]), "max": None},
            "y": {"min": None, "mid": float(df.loc[i, "y"]), "max": None},
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
            "sys": float(df.loc[i, "sys"]),
            "sys_beam": 3.7
            * data_central[idx]
            * 1e-2,  # 3.7% Normalization beam uncertainties
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
        "sys_beam": {
            "description": "systematic uncertainty due to beam Normalization",
            "treatment": "MULT",
            "type": "E143BEAM",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    fnames = glob.glob("rawdata/*.yaml")
    df = read_data(fnames)
    write_data(df)
