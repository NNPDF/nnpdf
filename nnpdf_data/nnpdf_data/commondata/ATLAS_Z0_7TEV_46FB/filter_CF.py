"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils_CF import get_data_values
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_7TEV_46FB_CF_data_central():
    """
    This function writes the central values to yaml files.
    """
    central_values = list(get_data_values())

    data_central_yaml = {"data_central": central_values}

    # write central values and kinematics to yaml file
    with open("data_CF.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_46FB_CF_data_central()
