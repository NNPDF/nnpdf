import yaml

Q2_POS = 5.0  # GeV^2

XVALUES_POS = [
    4.9999999999999998e-07,
    1.9407667236782136e-06,
    7.5331509514733370e-06,
    2.9240177382128657e-05,
    1.1349672651536727e-04,
    4.4054134013486355e-04,
    1.7099759466766963e-03,
    6.6373288312005724e-03,
    2.5763013859408150e-02,
    1.0000000000000001e-01,
    1.7999999999999999e-01,
    2.6000000000000001e-01,
    3.3999999999999997e-01,
    4.2000000000000004e-01,
    5.0000000000000000e-01,
    5.7999999999999996e-01,
    6.6000000000000003e-01,
    7.3999999999999999e-01,
    8.1999999999999995e-01,
    9.0000000000000002e-01,
]


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
