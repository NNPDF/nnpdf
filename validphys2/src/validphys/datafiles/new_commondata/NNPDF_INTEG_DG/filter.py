import yaml

Q2_INTEG = 1.0  # GeV^2

XVALUES_POS = [1e-5]


def write_data():
    kin = []
    for xvalue in XVALUES_POS:
        kin_value = {
            "x": {"min": None, "mid": xvalue, "max": None},
            "Q2": {"min": None, "mid": Q2_INTEG, "max": None},
            "y": {"min": None, "mid": 0.0, "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


if __name__ == "__main__":
    write_data()
