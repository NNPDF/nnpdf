"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_7TEV_46FB_data_central():
    """
    This function writes the central values to yaml files.
    """

    central_values_cc, central_values_cf = get_data_values()

    data_central_yaml_cc = {"data_central": central_values_cc}
    data_central_yaml_cf = {"data_central": central_values_cf}

    # write central values and kinematics to yaml file
    with open("data_cc.yaml", "w") as file:
        yaml.dump(data_central_yaml_cc, file, sort_keys=False)

    with open("data_cf.yaml", "w") as file:
        yaml.dump(data_central_yaml_cf, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_46FB_data_central()
