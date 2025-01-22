"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values, get_kinematics, get_systematics

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_7TEV_36FB_data_kinematic():
    """
    This function writes the systematics to yaml files.
    """

    central_values = get_data_values()

    kin = get_kinematics()

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_7TEV_36FB_systematics():
    """
    This function writes the systematics to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics()

    # error definition
    error_definitions = {}
    errors = []
    counter_1 = 1
    counter_2 = 0
    for sys in systematics:
        if sys[0]['name'] == 'stat':
            error_definitions[sys[0]['name']] = {
                "description": "Uncorrelated statistical uncertainties",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif sys[0]['name'] == 'uncor':
            error_definitions[sys[0]['name']] = {
                "description": f"Sys uncertainty idx: {counter_1}",
                "treatment": "MULT",
                "type": "UNCORR",
            }
            counter_1 += 1

        elif sys[0]['name'] == 'atlaslumi10':
            error_definitions[sys[0]['name']] = {
                "description": f"Sys uncertainty idx: {counter_1}",
                "treatment": "MULT",
                "type": "ATLASLUMI10",
            }
            counter_1 += 1

        else:
            error_definitions[sys[0]['name']] = {
                "description": f"Sys uncertainty idx: {counter_1}",
                "treatment": "MULT",
                "type": f"ATLASWZRAP36PB_{counter_2}",
            }
            counter_1 += 1
            counter_2 += 1

    for i in range(metadata['implemented_observables'][0]['ndata']):
        error_value = {}

        for sys in systematics:
            error_value[sys[0]['name']] = float(sys[0]['values'][i])

        errors.append(error_value)

    uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

    # write uncertainties
    with open(f"uncertainties.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_36FB_data_kinematic()
    filter_ATLAS_Z0_7TEV_36FB_systematics()
