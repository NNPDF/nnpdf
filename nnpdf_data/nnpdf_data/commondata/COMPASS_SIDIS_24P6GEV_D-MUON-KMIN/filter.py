import yaml
import os

# m_hadron = m_D = 1.876 GeV
E_cm = 24.6  # sqrt(s)=sqrt(2p_lep*m_hadron+m_hadron^2)=sqrt(2*160*1.876+1.876^2)=  24.6 GeV


def read_data(folder_path: str, tables: list):
    """
    htype: 0 for proton, 1 for deutron
    """
    collected_data = dict()
    collected_data["values"] = list()
    collected_data["errors_stat"] = list()
    collected_data["errors_sys"] = list()
    collected_data["kinematics_x"] = list()
    collected_data["kinematics_y"] = list()
    collected_data["kinematics_z"] = list()
    collected_data["kinematics_Q2"] = list()

    metadata_dict = {"tables": tables, "ndata_points": list()}

    for table in tables:
        with open(folder_path + f"data{table}.yaml", "r", encoding="utf-8") as file:
            file_dict = yaml.safe_load(file)
        values = file_dict["dependent_variables"][0]["values"]
        n_values = len(values)
        metadata_dict["ndata_points"].append(n_values)
        for i in range(n_values):
            collected_data["values"] = collected_data["values"] + [values[i]["value"]]
            collected_data["errors_stat"] = collected_data["errors_stat"] + [
                values[i]["errors"][0]["symerror"]
            ]
            collected_data["errors_sys"] = collected_data["errors_sys"] + [
                values[i]["errors"][1]["symerror"]
            ]
            collected_data["kinematics_Q2"] = collected_data["kinematics_Q2"] + [
                file_dict["independent_variables"][2]["values"][i]["value"]
            ]
            collected_data["kinematics_z"] = collected_data["kinematics_z"] + [
                [
                    file_dict["independent_variables"][3]["values"][i]["low"],
                    file_dict["independent_variables"][3]["values"][i]["value"],
                    file_dict["independent_variables"][3]["values"][i]["high"],
                ]
            ]
            collected_data["kinematics_x"] = collected_data["kinematics_x"] + [
                [
                    file_dict["independent_variables"][0]["values"][i]["low"],
                    file_dict["independent_variables"][0]["values"][i]["value"],
                    file_dict["independent_variables"][0]["values"][i]["high"],
                ]
            ]
            collected_data["kinematics_y"] = collected_data["kinematics_y"] + [
                [
                    file_dict["independent_variables"][1]["values"][i]["low"],
                    file_dict["independent_variables"][1]["values"][i]["value"],
                    file_dict["independent_variables"][1]["values"][i]["high"],
                ]
            ]
    return collected_data, metadata_dict


def write_data(collected_data: dict, folder_path: str):
    data_central_yaml = {"data_central": collected_data["values"]}
    with open(folder_path + f"data.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    n_items = len(collected_data["values"])

    # Write kin file
    kin = []
    for i in range(n_items):
        kin_value = {
            "x": {
                "min": collected_data["kinematics_x"][i][0],
                "mid": collected_data["kinematics_x"][i][1],
                "max": collected_data["kinematics_x"][i][2],
            },
            "y": {
                "min": collected_data["kinematics_y"][i][0],
                "mid": collected_data["kinematics_y"][i][1],
                "max": collected_data["kinematics_y"][i][2],
            },
            "z": {
                "min": collected_data["kinematics_z"][i][0],
                "mid": collected_data["kinematics_z"][i][1],
                "max": collected_data["kinematics_z"][i][2],
            },
            "Q2": {"min": None, "mid": collected_data["kinematics_Q2"][i], "max": None},
            "sqrts": {"min": None, "mid": E_cm, "max": None},
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}

    with open(folder_path + f"kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for i in range(n_items):
        # here uncertainties are symmetric
        e = {
            "stat": collected_data["errors_stat"][i],
            "sys": collected_data["errors_sys"][i],
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

    with open(folder_path + f"uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":

    # relative folder path
    folder_path = os.path.dirname(os.path.abspath(__file__)) + "/"

    # TODO Create a dict for easy running
    naming_dict = {
        # "KaPHadron-MLTP": [1],
        "KaMHadron-MLTP": [2],
    }

    # Wp
    for name, tables in naming_dict.items():
        collected_data, metadata_dict = read_data(folder_path + "rawdata/", tables)
        print(name.split("-")[0].lower(), metadata_dict)
        write_data(
            collected_data,
            folder_path=folder_path,
        )
