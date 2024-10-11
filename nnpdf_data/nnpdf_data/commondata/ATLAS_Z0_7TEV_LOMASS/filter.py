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


if __name__ == "__main__":
    filter_ATLAS_Z0_7TEV_LOMASS_data_kinetic()