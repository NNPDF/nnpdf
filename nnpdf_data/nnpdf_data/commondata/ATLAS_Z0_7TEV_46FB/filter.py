"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values, get_kinematics, get_systematics

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_7TEV_46FB_data_kinematic():
    """
    This function writes the systematics to yaml files.
    """

    central_values_cc, central_values_cf = get_data_values()

    kin_cc, kin_cf = get_kinematics()

    data_central_yaml_cc = {"data_central": central_values_cc}
    data_central_yaml_cf = {"data_central": central_values_cf}

    kinematics_yaml_cc = {"bins": kin_cc}
    kinematics_yaml_cf = {"bins": kin_cf}

    # write central values and kinematics to yaml file
    with open("data_cc.yaml", "w") as file:
        yaml.dump(data_central_yaml_cc, file, sort_keys=False)

    with open("data_cf.yaml", "w") as file:
        yaml.dump(data_central_yaml_cf, file, sort_keys=False)

    with open("kinematics_cc.yaml", "w") as file:
        yaml.dump(kinematics_yaml_cc, file, sort_keys=False)

    with open("kinematics_cf.yaml", "w") as file:
        yaml.dump(kinematics_yaml_cf, file, sort_keys=False)


def filter_ATLAS_Z0_7TEV_46FB_systematics():
    """
    This function writes the systematics to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics_cc, systematics_cf = get_systematics()

    # error definition
    error_definitions_cc = {}
    errors_cc = []

    error_definitions_cf = {}
    errors_cf = []

    counter = 1
    counter_3 = 0
    for sys in systematics_cc:

        if sys[0]['name'] == 'stat':
            error_definitions_cc[sys[0]['name']] = {
                "description": "Uncorrelated statistical uncertainties",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif sys[0]['name'] == 'AtlasLumi2011':
            error_definitions_cc[sys[0]['name']] = {
                "description": "'Sys uncertainty idx: 132'",
                "treatment": "MULT",
                "type": "ATLASLUMI11",
            }

        elif sys[0]['name'] == 'uncor' or sys[0]['name'] == 'uncor.1':
            counter_3 += 1

            if counter_3 == 1:
                counter += 1
                error_definitions_cc['sys_corr_1'] = {
                    "description": "Sys uncertainty idx: 1",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }
            elif counter_3 == 2:

                error_definitions_cc[sys[0]['name']] = {
                    "description": "Sys uncertainty idx: 133",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }

        else:
            error_definitions_cc['sys_corr_' + str(counter)] = {
                "description": "Sys uncertainty idx: " + str(counter),
                "treatment": "MULT",
                "type": f"{sys[0]['name']}",
            }
            counter += 1

    for i in range(metadata['implemented_observables'][0]['ndata']):
        error_value_cc = {}
        counter_2 = 0
        for sys in systematics_cc:
            if counter_2 == 0:
                error_value_cc[sys[0]['name']] = float(sys[0]['values'][i])
            else:
                error_value_cc['sys_corr_' + str(counter_2)] = float(sys[0]['values'][i])
            counter_2 += 1

        errors_cc.append(error_value_cc)

    uncertainties_yaml_cc = {"definitions": error_definitions_cc, "bins": errors_cc}

    # write uncertainties
    with open(f"uncertainties_cc.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml_cc, file, sort_keys=False)

    counter = 1
    counter_3 = 0
    for sys in systematics_cf:

        if sys[0]['name'] == 'stat':
            error_definitions_cf[sys[0]['name']] = {
                "description": "Uncorrelated statistical uncertainties",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif sys[0]['name'] == 'AtlasLumi2011':
            error_definitions_cf[sys[0]['name']] = {
                "description": "'Sys uncertainty idx: 132'",
                "treatment": "MULT",
                "type": "ATLASLUMI11",
            }

        elif sys[0]['name'] == 'uncor' or sys[0]['name'] == 'uncor.1':
            counter_3 += 1

            if counter_3 == 1:
                counter += 1
                error_definitions_cf['sys_corr_1'] = {
                    "description": "Sys uncertainty idx: 1",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }
            elif counter_3 == 2:

                error_definitions_cf[sys[0]['name']] = {
                    "description": "Sys uncertainty idx: 133",
                    "treatment": "MULT",
                    "type": "UNCORR",
                }

        else:
            error_definitions_cf['sys_corr_' + str(counter)] = {
                "description": "Sys uncertainty idx: " + str(counter),
                "treatment": "MULT",
                "type": f"{sys[0]['name']}",
            }
            counter += 1

    for i in range(metadata['implemented_observables'][1]['ndata']):
        error_value_cf = {}
        counter_2 = 0
        for sys in systematics_cf:
            if counter_2 == 0:
                error_value_cf[sys[0]['name']] = float(sys[0]['values'][i])
            else:
                error_value_cf['sys_corr_' + str(counter_2)] = float(sys[0]['values'][i])
            counter_2 += 1

        errors_cf.append(error_value_cf)

    uncertainties_yaml_cf = {"definitions": error_definitions_cf, "bins": errors_cf}

    # write uncertainties
    with open(f"uncertainties_cf.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml_cf, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_46FB_data_kinematic()
    filter_ATLAS_Z0_7TEV_46FB_systematics()
