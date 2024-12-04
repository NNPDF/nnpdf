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


if __name__ == "__main__":
    get_data_values()
