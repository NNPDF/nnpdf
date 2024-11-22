import yaml
import os

ECM = 10.3  # sqrt(s)=sqrt(2p_lep*m_hadron+m_hadron^2)=sqrt(2*27.6*1.876+1.876^2)= 10.3 in GeV


def read_data(folder_path: str, tables: list, htype: int):
    """
    htype: 0 for proton, 1 for deutron
    """
    collected_data = dict()
    collected_data["values"] = list()
    collected_data["errors_stat"] = list()
    collected_data["errors_sys"] = list()
    collected_data["kinematics_x"] = list()
    collected_data["kinematics_z"] = list()

    metadata_dict = {"htype": htype, "tables": tables, "ndata_points": list()}

    for table in tables:
        with open(folder_path + f"Table{table}.yaml", "r", encoding="utf-8") as file:
            file_dict = yaml.safe_load(file)
        z_str = file_dict["dependent_variables"][htype]["qualifiers"][2]["value"]
        z_min, z_max = map(float, z_str.split("-"))
        values = file_dict["dependent_variables"][htype]["values"]
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
            collected_data["kinematics_x"] = collected_data["kinematics_x"] + [
                file_dict["independent_variables"][htype]["values"][i]["value"]
            ]
            collected_data["kinematics_z"] = collected_data["kinematics_z"] + [
                [
                    z_min,
                    z_max,
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
            "z": {
                "min": collected_data["kinematics_z"][i][0],
                "mid": None,
                "max": collected_data["kinematics_z"][i][1],
            },
            "x": {"min": None, "mid": collected_data["kinematics_x"][i], "max": None},
            "sqrts": {"min": None, "mid": ECM, "max": None},
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

    # Get the path of the current file
    folder_path = os.path.dirname(os.path.abspath(__file__)) + "/"

    # TODO Create a dict for easy running
    naming_dict = {
        # "PiPProton-MLTP": [33, 34, 35, 36],
        # "PiPDeutron-MLTP": [33, 34, 35, 36],
        # "PiMProton-MLTP": [37, 38, 39, 40],
        "PiMDeutron-MLTP": [37, 38, 39, 40],
        # "KaMProton-MLTP": [45, 46, 47, 48],
        # "KaMDeutron-MLTP": [45, 46, 47, 48],
        # "KaPProton-MLTP": [41, 42, 43, 44],
        # "KaPDeutron-MLTP": [41, 42, 43, 44],
    }

    # Wp
    for name, tables in naming_dict.items():
        if "Proton" in name:
            htype = 0
        else:
            htype = 1

        if name.upper() in folder_path:
            a = 1

        collected_data, metadata_dict = read_data(
            folder_path + "rawdata/", tables, htype
        )
        print(name.split("-")[0].lower(), metadata_dict)
        write_data(
            collected_data,
            folder_path=folder_path,
        )
