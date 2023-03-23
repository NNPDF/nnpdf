# Filter for DYE866R

import yaml


def filter_DYE866P():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]

    data_central = []
    kin = []
    error = []

    for i in tables:
        hepdata_tables = (
            "rawdata/HEPData-ins613362-v" + str(version) + "-Table_" + str(i) + ".yaml"
        )
        with open(hepdata_tables, "r") as file:
            input = yaml.safe_load(file)
        values = input["dependent_variables"][0]["values"]
        sqrts = float(input["dependent_variables"][0]["qualifiers"][2]["value"])

        for j in range(len(values)):

            data_central_value = values[j]["value"]
            data_central.append(data_central_value)
            xF = input["independent_variables"][1]["values"][j]["value"]
            M_high = input["independent_variables"][0]["values"][j]["high"]
            M_low = input["independent_variables"][0]["values"][j]["low"]
            # M = input["independent_variables"][0]["values"][j]["value"]
            # value for central M is missing for one point, so I m computing it.
            # This is the way it is done in old implementation
            M = 0.5 * (M_high + M_low)

            kin_value = {
                "xF": {"min": None, "mid": xF, "max": None},
                "M": {"min": M_low, "mid": M, "max": M_high},
                "sqrts": {"min": None, "mid": sqrts, "max": None},
            }
            kin.append(kin_value)

            error_value = {
                "stat_1": input["dependent_variables"][0]["values"][j]["errors"][0][
                    "symerror"
                ],
                "syst_1": input["dependent_variables"][0]["values"][j]["errors"][1][
                    "symerror"
                ],
                "syst_2": data_central_value * 6.5 * 1e-2,
            }
            error.append(error_value)

    error_definition = {
        "stat_1": {
            "description": "total statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "syst_1": {
            "description": "total systematic uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "syst_2": {
            "description": "Additional systematic uncertainty in the normalization",
            "treatment": "MULT",
            "type": "CORR",
        },
    }

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}
    uncertainties_yaml = {"definition": error_definition, "bins": error}

    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


filter_DYE866P()
