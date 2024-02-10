import yaml

Q2_POS = 1.0  # GeV^2

XVALUES_POS = [
    0.001,
    0.0018059880241530946,
    0.003251185184362275,
    0.005820073691256094,
    0.010318697769338017,
    0.018005236582438074,
    0.03064850962450646,
    0.050358677992155945,
    0.07909372187858188,
    0.11804121183773833,
    0.16730811062687573,
    0.2261040353427545,
    0.2931760066717947,
    0.3671850159384811,
    0.44691296552565674,
    0.5313348186291872,
    0.6196196728514175,
    0.7111045143674228,
    0.8052623301561946,
    0.9016730963186919,
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
