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

def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys

def get_uncertainties():
    """
    returns the uncertainties.
    """

    ndat = 5
    # Produce covmat of form [[W-/W+],[0],
    #                        [0],[W-*/W+*]]
    covmat = np.zeros((4*ndat, 4*ndat)) # Multiply by 4 because of W+/- and */not *
   
    def edit_covmat(filename, offset):
        with open(filename) as f:
            data = yaml.safe_load(f)
        flat_values = [v["value"] for v in data["dependent_variables"][0]["values"]]
        matrix = np.array(flat_values).reshape((2 * ndat, 2 * ndat))
        covmat[offset:offset + 2 * ndat, offset:offset + 2 * ndat] = matrix

    edit_covmat("rawdata/HEPData-ins2628732-v1-Table_16.yaml", offset=0)
    edit_covmat("rawdata/HEPData-ins2628732-v1-Table_18.yaml", offset=2 * ndat)

    sys = decompose_covmat(covmat)
