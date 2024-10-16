"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml


def get_data_values():
    """
    returns the central data values in the form of a list.
    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins1234228-v1-Table_1.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        # store data central and convert the units
        data_central.append(value['value'] * 1000)

    return data_central


if __name__ == "__main__":
    get_data_values()
