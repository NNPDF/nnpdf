import yaml
from filter_utils import (
                            get_data_values, get_kinematics
                        )

def filter_ATLAS_1JET_8TEV_data_kinetic():
    """
    write kinematics and central data values
    in kinematics.yaml and data.yaml files
    respectively.
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    # get kinematics from hepdata tables
    kin = get_kinematics(tables,version)

    # get central values from hepdata tables
    data_central = get_data_values(tables,version)
    
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


if __name__ == "__main__" :

    # write kinematics and central data values
    filter_ATLAS_1JET_8TEV_data_kinetic()