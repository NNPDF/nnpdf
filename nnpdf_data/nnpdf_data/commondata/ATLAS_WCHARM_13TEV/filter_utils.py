"""
This module contains helper functions that are used to extract the data values 
from the rawdata files.
"""

import yaml
import pandas as pd
import numpy as np
from numpy.linalg import eig


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

def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):
    r"""Convert the covariance matrix to a matrix of
    artificial uncertainties.

    NOTE: This function has been taken from validphys.newcommondata_utils.
    If those utils get merged in the future, we can replace this.

    Parameters
    ----------
    ndata : integer
        Number of data points
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.
    no_of_norm_mat : int
        Normalized covariance matrices may have an eigenvalue
        of 0 due to the last data point not being linearly
        independent. To allow for this, the user should input
        the number of normalized matrices that are being treated
        in an instance. For example, if a single covariance matrix
        of a normalized distribution is being processed, the input
        would be 1. If a covariance matrix contains pertains to
        3 normalized datasets (i.e. cross covmat for 3
        distributions), the input would be 3. The default value is
        0 for when the covariance matrix pertains to an absolute
        distribution.

    Returns
    -------
    artunc : list
        A two dimensional matrix (given as a list of lists)
        which contains artificial uncertainties to be added
        to the commondata. i^th row (or list) contains the
        artificial uncertainties of the i^th data point.

    """
    epsilon = -0.0000000001
    neg_eval_count = 0
    psd_check = True
    covmat = np.zeros((ndata, ndata))
    artunc = np.zeros((ndata, ndata))
    for i in range(len(covmat_list)):
        a = i // ndata
        b = i % ndata
        covmat[a][b] = covmat_list[i]
    eigval, eigvec = eig(covmat)
    for j in range(len(eigval)):
        if eigval[j] < epsilon:
            psd_check = False
        elif eigval[j] > epsilon and eigval[j] <= 0:
            neg_eval_count = neg_eval_count + 1
            if neg_eval_count == (no_of_norm_mat + 1):
                psd_check = False
        elif eigval[j] > 0:
            continue
    if psd_check == False:
        raise ValueError("The covariance matrix is not positive-semidefinite")
    else:
        for i in range(ndata):
            for j in range(ndata):
                if eigval[j] < 0:
                    continue
                else:
                    artunc[i][j] = eigvec[i][j] * np.sqrt(eigval[j])
    return artunc.tolist()

def get_artificial_uncertainties():
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

    covmat_list = covmat.flatten().tolist()
    artificial_sys = np.array(covmat_to_artunc(4*ndat, covmat_list))
    uncertainties = []
    uncertainties.append([{"name": "stat", "values": np.zeros(4*ndat)}])

    for i in range(len(artificial_sys)):
        name = f"sys_{i}"
        values = artificial_sys[:,i]
        uncertainties.append([{"name": name, "values": values}])
    
    return uncertainties
    return artificial_sys

