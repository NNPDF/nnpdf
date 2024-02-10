import yaml

Q2_POS = 1.0  # GeV^2

XVALUES_POS = [1e-5]


def write_data():
    data_central = []
    for _ in range(len(XVALUES_POS)):
        data_central.append(0.0)

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for xvalue in XVALUES_POS:
        kin_value = {
            "x": {"min": None, "mid": xvalue, "max": None},
            "Q2": {"min": None, "mid": Q2_POS, "max": None},
            "y": {"min": None, "mid": 0.0, "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Write unc file
    error = []
    for _ in range(len(XVALUES_POS)):
        e = {
            "stat": 0.0,
        }
        error.append(e)

    error_definition = {
        "stat": {
            "description": "statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
    }

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    write_data()
