"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values, get_kinematics, get_systematics
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_7TEV_46FB_CC_data_kinematics():
    """
    This function writes the central values to yaml files.
    """
    central_values = list(get_data_values())

    kin = get_kinematics()

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data_CC.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics_CC.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_7TEV_46FB_CC_systematics():
    """
    This function writes the systematics to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics()

    # error definition
    error_definitions = {}
    errors = []

    counter = 1
    counter_3 = 0
    for sys in systematics:
        print(sys[0]['name'])
        if sys[0]['name'] == 'stat':
            error_definitions[sys[0]['name']] = {
                "description": "Uncorrelated statistical uncertainties",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif sys[0]['name'] == 'AtlasLumi2011':
            error_definitions[sys[0]['name']] = {
                "description": "'Sys uncertainty idx: 132'",
                "treatment": "MULT",
                "type": "ATLASLUMI11",
            }

        elif sys[0]['name'] == 'uncor' or sys[0]['name'] == 'uncor.1':
            counter_3 += 1
            print(counter_3)
            if counter_3 == 1:
                counter += 1
                error_definitions['sys_corr_1'] = {
                    "description": "Sys uncertainty idx: 1",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }
            elif counter_3 == 2:
                print('yes',sys[0]['name'])
                error_definitions[sys[0]['name']] = {
                    "description": "Sys uncertainty idx: 133",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }

        else:
            error_definitions['sys_corr_' + str(counter)] = {
                "description": "Sys uncertainty idx: " + str(counter),
                "treatment": "MULT",
                "type": f"{sys[0]['name']}",
            }
            counter += 1

    for i in range(metadata['implemented_observables'][0]['ndata']):
        error_value = {}
        counter_2 = 0
        for sys in systematics:
            if counter_2 == 0:
                error_value[sys[0]['name']] = float(sys[0]['values'][i])
            else:
                error_value['sys_corr_' + str(counter_2)] = float(sys[0]['values'][i])
            counter_2 += 1

        errors.append(error_value)

    uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

    # write uncertainties
    with open(f"uncertainties.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_46FB_CC_data_kinematics()
    filter_ATLAS_Z0_7TEV_46FB_CC_systematics()
