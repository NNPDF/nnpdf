# Filter for NMCPD

import yaml


def filter_E605():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]
    npoints = metadata["hepdata"]["npoints"]

    data_central = []
    kin = []
    error = []
    for i, n in zip(tables, npoints):
        hepdata_tables = (
            "rawdata/HEPData-ins302822-v" + str(version) + "-Table_" + str(i) + ".yaml"
        )
        with open(hepdata_tables, "r") as file:
            input = yaml.safe_load(file)

        y = float(input["dependent_variables"][0]["qualifiers"][2]["value"])
        sqrts = float(input["dependent_variables"][0]["qualifiers"][1]["value"])

        for j in range(n):

            data_central_value = input["dependent_variables"][0]["values"][j]["value"]
            data_central.append(data_central_value)
            sqrtau = input["independent_variables"][0]["values"][j]["value"]
            kin_value = {
                "s": {"min": None, "mid": sqrts**2, "max": None},
                "sqrtau": {"min": None, "mid": sqrtau, "max": None},
                "y": {"min": None, "mid": y, "max": None},
            }
            kin.append(kin_value)
            error_value = {
                "stat_1": input["dependent_variables"][0]["values"][j]["errors"][0][
                    "symerror"
                ],
                "syst_1": float(
                    input["dependent_variables"][0]["values"][j]["errors"][1][
                        "symerror"
                    ].rstrip("%")
                )
                * data_central_value
                * 1e-2,
                "syst_2": float(
                    input["dependent_variables"][0]["values"][j]["errors"][2][
                        "symerror"
                    ].rstrip("%")
                )
                * data_central_value
                * 1e-2,
            }
            error.append(error_value)

    error_definition = {
        "stat_1": {
            "description": "total statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "syst_1": {
            "description": "normalization systematic uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
        "syst_2": {
            "description": "total systematic uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
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


filter_E605()
