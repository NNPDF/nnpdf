import yaml
import numpy as np


def filter_DYE866R():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    obs = metadata["implemented_observables"][0]
    tables = obs["tables"]

    data_central = []
    kin = []
    error = []

    for i in tables:
        hepdata_tables = (
            "rawdata/HEPData-ins554316-v" + str(version) + "-Table_" + str(i) + ".yaml"
        )
        with open(hepdata_tables, "r") as file:
            input = yaml.safe_load(file)

        values = input["dependent_variables"][0]["values"]
        sqrts = float(input["dependent_variables"][0]["qualifiers"][3]["value"])

        for j in range(len(values)):

            data_central_value = values[j]["value"]
            data_central.append(data_central_value)
            xF = input["independent_variables"][1]["values"][j]["value"]
            M = input["independent_variables"][3]["values"][j]["value"]
            J = np.sqrt(xF**2 + 4 * M**2 / sqrts**2)
            y = float(0.5 * np.log((J + xF) / (J - xF)))

            kin_value = {
                "y": {"min": None, "mid": y, "max": None},
                "m2": {"min": None, "mid": M**2, "max": None},
                "sqrts": {"min": None, "mid": sqrts, "max": None},
            }
            kin.append(kin_value)

            error_value = {
                "stat_1": input["dependent_variables"][0]["values"][j]["errors"][0][
                    "symerror"
                ],
                "syst_1": data_central_value * 0.97 * 1e-2,
            }
            error.append(error_value)

    error_definition = {
        "stat_1": {
            "description": "total statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "syst_1": {
            "description": "total systematic uncertainty, 0.97%",
            "treatment": "ADD",
            "type": "CORR",
        },
    }

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}
    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


filter_DYE866R()
