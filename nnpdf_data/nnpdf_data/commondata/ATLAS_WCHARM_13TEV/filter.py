"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import (
    get_data_values,
    get_kinematics,
    get_artificial_uncertainties,
    get_uncertainties,
)
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_WCHARM_13TEV_data_kinematic():
    """
    This function writes the systematics to yaml files.
    """

    central_values = get_data_values()

    kin = get_kinematics()

    data_central_yaml_WMDP = {"data_central": central_values[:5]}
    data_central_yaml_WPDM = {"data_central": central_values[5:10]}
    data_central_yaml_WMDP_star = {"data_central": central_values[10:15]}
    data_central_yaml_WPDM_star = {"data_central": central_values[15:20]}

    kinematics_yaml_WMDP = {"bins": kin[:5]}
    kinematics_yaml_WPDM = {"bins": kin[5:10]}
    kinematics_yaml_WMDP_star = {"bins": kin[10:15]}
    kinematics_yaml_WPDM_star = {"bins": kin[15:20]}

    # write central values and kinematics to yaml file
    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        with open(f"data_{channel}.yaml", "w") as file:
            data_central_yaml = eval(f"data_central_yaml_{channel}")
            yaml.dump(data_central_yaml, file, sort_keys=False)

    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        kinematics_yaml = eval(f"kinematics_yaml_{channel}")
        with open(f"kinematics_{channel}.yaml", "w") as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_get_artificial_uncertainties():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics_full = get_artificial_uncertainties()

    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:

        systematics = systematics_full[channel]

        error_definitions = {}
        errors = []

        for sys in systematics:
            if sys['name'] == 'stat':
                error_definitions[sys['name']] = {
                    "description": "Uncorrelated statistical uncertainties",
                    "treatment": "ADD",
                    "type": "UNCORR",
                }
            else:
                error_definitions[sys['name']] = {
                    "description": "Systematic uncertainty",
                    "treatment": "MULT",  # Not sure if this is correct
                    "type": "CORR",
                }

        for i in range(metadata['implemented_observables'][0]['ndata']):
            error_value = {}

            for sys in systematics:
                error_value[sys['name']] = float(sys['values'][i])

            errors.append(error_value)

        uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

        # write uncertainties
        with open(f"uncertainties_covariances_{channel}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


def filter_get_systematics():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics_full, deltas = get_uncertainties()
    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        systematics = systematics_full[channel]
        error_definitions = {}
        errors = []
        # print("Systematics:", systematics)
        for sys in systematics:
            if sys['name'] == 'stat':
                error_definitions[sys['name']] = {
                    "description": "Uncorrelated statistical uncertainties",
                    "treatment": "ADD",
                    "type": "UNCORR",
                }
            else:
                error_definitions[sys['name']] = {
                    "description": "Systematic uncertainty",
                    "treatment": "MULT",  # Not sure if this is correct
                    "type": "CORR",
                }

        for i in range(metadata['implemented_observables'][0]['ndata']):
            error_value = {}

            for sys in systematics:

                error_value[sys['name']] = float(sys['values'][i])

            errors.append(error_value)

        uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

        # write uncertainties
        with open(f"uncertainties_{channel}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_WCHARM_13TEV_data_kinematic()
    filter_get_artificial_uncertainties()
    filter_get_systematics()
