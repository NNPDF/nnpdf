"""
filter.py module for ATLAS_WPWM_13TEV dataset
When running `python filter.py` the relevant uncertainties , data and kinematics yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_13TEV` directory.
"""

import yaml
from filter_utils import get_data_values, get_systematics
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

MZ2 = 91.1876**2


def filter_ATLAS_Z0_13TEV_TOT_data_kinetic():
    """
    This function writes the central values and kinematics to yaml files.
    """

    kin = [
        {
            'm_Z2': {'min': None, 'mid': MZ2, 'max': None},
            'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
        }
    ]

    # only keep the last entry corresponding to Z observable
    central_values = [get_data_values()[-1]]

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_13TEV_TOT_systematics():
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
        if sys[0]['name'] == 'stat':
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif sys[0]['name'] == 'ATLAS_LUMI':
            error_definitions["ATLASLUMI15PART"] = {
                "description": f"ATLASLUMI15PART",
                "treatment": "MULT",
                "type": "ATLASLUMI15PART",
            }

        else:
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "ADD",
                "type": f"{sys[0]['name']}",
            }

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
    filter_ATLAS_Z0_13TEV_TOT_data_kinetic()
    filter_ATLAS_Z0_13TEV_TOT_systematics()
