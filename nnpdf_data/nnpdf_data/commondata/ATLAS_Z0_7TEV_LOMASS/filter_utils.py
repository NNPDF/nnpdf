"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml
import pandas as pd
import numpy as np


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins1288706-v1-Table_6.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    y_values = [12.41, 22.57, 14.64, 6.73, 2.81, 1.27]

    for i, M in enumerate(input["independent_variables"][0]['values']):

        kin_value = {
            'y': {'min': None, 'mid': y_values[i], 'max': None},
            'M2': {'min': M['low'] ** 2, 'mid': 0.5 * (M['low'] + M['high']), 'max': M['high']},
            'sqrts': {'min': None, 'mid': 7000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values():
    """
    returns the central data values in the form of a numpy array.
    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins1288706-v1-Table_6.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        # store data central and convert the units
        data_central.append(value['value'] * 1000)

    return np.array(data_central)


def get_systematics():
    """
    returns the absolute systematic uncertainties in the form of a pandas dataframe.
    """
    sys_rawdata_path = "rawdata/ATLASLOMASSDY11EXT.csv"

    df = pd.read_csv(sys_rawdata_path)
    data_central = get_data_values()

    # convert (MULT) percentage unc to absolute unc
    abs_unc_df = (df.T[3:] * data_central).T / 100

    return abs_unc_df
