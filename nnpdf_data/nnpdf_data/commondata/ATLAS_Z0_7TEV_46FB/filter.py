"""
When running `python filter.py` the relevant data yaml
file will be created in the `nnpdf_data/commondata/ATLAS_WPWM_7TEV_46FB` directory.
"""

import yaml
from filter_utils import get_data_values, get_kinematics

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


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_46FB_data_kinematic()
