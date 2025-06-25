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

    for i in range(19, 23):
        hepdata_table = f"rawdata/HEPData-ins2628732-v1-Table_{i}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']

        for value in values:
            # store data central and convert the units and apply the correction factor
            data_central.append(value['value'])

    return data_central


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    for i in range(19, 23):
        hepdata_table = f"rawdata/HEPData-ins2628732-v1-Table_{i}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        for i, M in enumerate(input["independent_variables"][0]['values']):
            kin_value = {
                'abs_eta': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
                'm_W2': {'min': None, 'mid': 6.46046213e03, 'max': None},
                'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
            }
            kin.append(kin_value)

    return kin
