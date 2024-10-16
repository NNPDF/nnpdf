"""
filter.py module for ATLAS_Z0_7TEV_49FB dataset
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_Z0_7TEV_LOMASS` directory.
"""

import yaml
from filter_utils import get_data_values


def filter_ATLAS_Z0_7TEV_49FB_data_central():
    """
    This function writes the central values to yaml files.
    """
    central_values = list(get_data_values())

    data_central_yaml = {"data_central": central_values}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_49FB_data_central()
