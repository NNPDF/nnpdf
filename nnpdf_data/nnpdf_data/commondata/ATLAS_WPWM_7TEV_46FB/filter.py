"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values, get_kinematics


def filter_ATLAS_WPWM_7TEV_46FB_data_kinematic():
    """
    This function writes the central values to yaml files.
    """
    central_values = list(get_data_values())

    kin = get_kinematics()

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_WPWM_7TEV_46FB_data_kinematic()