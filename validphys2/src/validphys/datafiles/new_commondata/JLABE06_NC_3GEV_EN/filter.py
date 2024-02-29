import pandas as pd
import yaml
import glob
from io import StringIO


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            # Get the data as string
            data = file.read()

            # Split the data into sections based on empty lines
            sections = data.strip().split('\n\n')

            # Get the last section, excluding the lines containing dashes and the two lines before them
            last_section = sections[-1].split('\n')[3:]

            # Correctly set the columns
            last_section[0] = "x A_1 A_stat A_sys G G_stat G_sys"

            # Convert the modified section into a Pandas DataFrame
            df_temp = pd.read_csv(StringIO('\n'.join(last_section)), delim_whitespace=True)

            # Set the Q^2 and y value
            df_temp["Q2"] = "3.08"
            df_temp["y"] = "0.0"

            # Concat with the total dataframe
            df = pd.concat([df, df_temp])

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
            "x": {"min": None, "mid": float(df.loc[i, "x"]), "max": None},
            "Q2": {"min": None, "mid": float(df.loc[i, "Q2"]), "max": None},
            "y": {"min": None, "mid": float(df.loc[i, "y"]), "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(len(df)):
        e = {
            "stat": float(df.loc[i, "G_stat"]),
            "sys": float(df.loc[i, "G_sys"]),
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

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    fnames = glob.glob("rawdata/*.txt")
    nnames = sorted([i for i in fnames])
    df = read_data(nnames)
    write_data(df)
