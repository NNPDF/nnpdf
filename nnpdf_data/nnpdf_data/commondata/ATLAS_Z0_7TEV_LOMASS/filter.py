"""
filter.py module for ATLAS_Z0_7TEV_LOMASS dataset

When running `python filter.py` the relevant uncertainties , data and kinematics yaml
file will be created in the `nnpdf_data/commondata/ATLAS_Z0_7TEV_LOMASS` directory.
"""

import yaml
from filter_utils import get_kinematics, get_data_values, get_systematics


def filter_ATLAS_Z0_7TEV_LOMASS_data_kinetic():
    """
    This function writes the central values and kinematics to yaml files.
    """

    kin = get_kinematics()
    central_values = list(get_data_values())

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_7TEV_LOMASS_systematics():
    """
    This function writes the systematics to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics()

    # error definition
    error_definitions = {}
    errors = []

    for sys in systematics:

        if (
            (sys[0]['name'] == 'stat')
            or (sys[0]['name'] == 'sys_res')
            or (sys[0]['name'] == 'sys_MC')
        ):
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "MULT",
                "type": "UNCORR",
            }

        else:
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "MULT",
                "type": "CORR",
            }

    #
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
    filter_ATLAS_Z0_7TEV_LOMASS_data_kinetic()
    filter_ATLAS_Z0_7TEV_LOMASS_systematics()
