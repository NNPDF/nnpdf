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

    tables = [5, 3]

    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins928289-v1-Table_{table}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']

        for value in values:
            # store data central and convert the units and apply the correction factor
            data_central.append(value['value'] * 1000 * 1.0187)

    return data_central


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    tables = [5, 3]

    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins928289-v1-Table_{table}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        for i, M in enumerate(input["independent_variables"][0]['values']):
            kin_value = {
                'abs_eta': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
                'm_W2': {'min': None, 'mid': 6463.838404, 'max': None},
                'sqrts': {'min': None, 'mid': 7000.0, 'max': None},
            }
            kin.append(kin_value)

    return kin


if __name__ == "__main__":
    get_data_values()
    get_kinematics()
