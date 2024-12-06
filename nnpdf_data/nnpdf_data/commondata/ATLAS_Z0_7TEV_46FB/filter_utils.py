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

    data_central_cc = []
    data_central_cf = []
    cc_tables = [11, 12, 13]
    cf_tables = [14, 15]

    for table in cc_tables:
        hepdata_table = f"rawdata/HEPData-ins1502620-v1-Table_{table}.yaml"
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for value in values:
            # store data central and convert the units
            data_central_cc.append(value['value'] * 1000)

    for table in cf_tables:
        hepdata_table = f"rawdata/HEPData-ins1502620-v1-Table_{table}.yaml"
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for value in values:
            # store data central and convert the units
            data_central_cf.append(value['value'] * 1000)

    return data_central_cc, data_central_cf


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin_cc = []
    kin_cf = []
    cc_tables = [11, 12, 13]
    cf_tables = [14, 15]

    # Define a mapping for table numbers to av_m_ll2 values
    av_m_ll2_mapping = {11: 56**2, 12: 91**2, 13: 133**2, 14: 91**2, 15: 133**2}

    for table in cc_tables:
        hepdata_table = f"rawdata/HEPData-ins1502620-v1-Table_{table}.yaml"
        av_m_ll2 = av_m_ll2_mapping[table]
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        for i, M in enumerate(input["independent_variables"][0]['values']):
            kin_value = {
                'abs_eta': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
                'm_ll2': {'min': None, 'mid': av_m_ll2, 'max': None},
                'sqrts': {'min': None, 'mid': 7000.0, 'max': None},
            }
            kin_cc.append(kin_value)

    for table in cf_tables:
        hepdata_table = f"rawdata/HEPData-ins1502620-v1-Table_{table}.yaml"
        av_m_ll2 = av_m_ll2_mapping[table]
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        for i, M in enumerate(input["independent_variables"][0]['values']):
            kin_value = {
                'abs_eta': {'min': None, 'mid': (0.5 * (M['low'] + M['high'])), 'max': None},
                'm_ll2': {'min': None, 'mid': av_m_ll2, 'max': None},
                'sqrts': {'min': None, 'mid': 7000.0, 'max': None},
            }
            kin_cf.append(kin_value)

    return kin_cc, kin_cf


if __name__ == "__main__":
    get_data_values()
    get_kinematics()
