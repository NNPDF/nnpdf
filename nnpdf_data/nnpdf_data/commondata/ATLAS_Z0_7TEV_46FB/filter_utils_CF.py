"""
This module contains helper functions that are used to extract the data values 
from the rawdata files.
"""

import yaml

import pandas as pd
import numpy as np


def get_data_values():
    """
    returns the central data values in the form of a list.
    """

    data_central = []

    hepdata_table_1 = f"rawdata/HEPData-ins1502620-v1-Table_14.yaml"
    hepdata_table_2 = f"rawdata/HEPData-ins1502620-v1-Table_15.yaml"

    with open(hepdata_table_1, 'r') as file:
        input_1 = yaml.safe_load(file)

    with open(hepdata_table_2, 'r') as file:
        input_2 = yaml.safe_load(file)

    values_1 = input_1['dependent_variables'][0]['values']
    values_2 = input_2['dependent_variables'][0]['values']

    values = values_1 + values_2

    for value in values:
        # store data central and convert the units
        data_central.append(value['value'] * 1000)

    return data_central


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    hepdata_table_1 = f"rawdata/HEPData-ins1502620-v1-Table_14.yaml"
    hepdata_table_2 = f"rawdata/HEPData-ins1502620-v1-Table_15.yaml"

    with open(hepdata_table_1, 'r') as file:
        input_1 = yaml.safe_load(file)

    with open(hepdata_table_2, 'r') as file:
        input_2 = yaml.safe_load(file)

    for i, M in enumerate(input_1["independent_variables"][0]['values']):

        kin_value = {
            'k1': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
            'k2': {'min': None, 'mid': 8281.0, 'max': None},
            'k3': {'min': None, 'mid': 7000.0, 'max': None},
        }

        kin.append(kin_value)

    for i, M in enumerate(input_2["independent_variables"][0]['values']):

        kin_value = {
            'k1': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
            'k2': {'min': None, 'mid': 17689.0, 'max': None},
            'k3': {'min': None, 'mid': 7000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


if __name__ == "__main__":
    get_data_values()
    get_kinematics()
