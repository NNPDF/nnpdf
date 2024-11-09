"""
This module contains helper functions that are used to extract the data values 
from the rawdata files.
"""

import yaml
# import pandas as pd
# import numpy as np


def get_data_values():
    """
    returns the central data values in the form of a list.
    """

    data_central = []

    hepdata_table_1 = f"rawdata/HEPData-ins1502620-v1-Table_11.yaml"
    hepdata_table_2 = f"rawdata/HEPData-ins1502620-v1-Table_12.yaml"
    hepdata_table_3 = f"rawdata/HEPData-ins1502620-v1-Table_13.yaml"

    with open(hepdata_table_1, 'r') as file:
        input_1 = yaml.safe_load(file)

    with open(hepdata_table_2, 'r') as file:
        input_2 = yaml.safe_load(file)

    with open(hepdata_table_3, 'r') as file:
        input_3 = yaml.safe_load(file)

    values_1 = input_1['dependent_variables'][0]['values']
    values_2 = input_2['dependent_variables'][0]['values']
    values_3 = input_3['dependent_variables'][0]['values']

    values = values_1 + values_2 + values_3

    for value in values:
        # store data central and convert the units
        data_central.append(value['value'] * 1000)

    return data_central


if __name__ == "__main__":
    get_data_values()
   