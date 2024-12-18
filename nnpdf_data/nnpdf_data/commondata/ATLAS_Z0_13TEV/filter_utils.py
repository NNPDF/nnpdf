"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml
import numpy as np
from nnpdf_data.filter_utils.utils import decompose_covmat


UNIT_CONVERSION = 1000000
TABLES = [9, 8, 11]  # order is W-, W+, Z


def get_data_values():
    """
    returns the central data values in the form of a list.
    """
    name_data = lambda tab: f"rawdata/HEPData-ins1436497-v1-Table_{tab}.yaml"

    data_central = []

    for tab in TABLES:
        with open(name_data(tab), 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        data_central.append(values[0]['value'] * UNIT_CONVERSION)

    return data_central


def get_uncertainties():
    """
    Returns array of shape (3,3)
    Each row corresponds to a different observable: (W-, W+, Z)
    Each column corresponds to a different systematic: (stat, sys, lumi)

    See table 3 of paper: https://arxiv.org/abs/1603.09222
    """

    name_data = lambda tab: f"rawdata/HEPData-ins1436497-v1-Table_{tab}.yaml"

    uncertainties = []

    for tab in TABLES:
        with open(name_data(tab), 'r') as file:
            input = yaml.safe_load(file)
        errors = input['dependent_variables'][0]['values'][0]['errors']
        uncertainties.append(
            np.array([errors[0]['symerror'], errors[1]['symerror'], errors[2]['symerror']])
        )

    return np.array(uncertainties) * UNIT_CONVERSION


def get_correlation_matrix():
    """
    See extra material page: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/STDM-2015-03/tabaux_03.pdf

    Note that this does not include the normalisation uncertainty due to the luminosity.
    """

    correlation_matrix = np.ones((3, 3))
    correlation_matrix[0, 1] = 0.93
    correlation_matrix[1, 0] = correlation_matrix[0, 1]
    correlation_matrix[0, 2] = 0.18
    correlation_matrix[2, 0] = correlation_matrix[0, 2]
    correlation_matrix[1, 2] = 0.19
    correlation_matrix[2, 1] = correlation_matrix[1, 2]

    return correlation_matrix


def get_covariance_matrices():
    """
    For the systematics see Table 3 of paper: https://arxiv.org/abs/1603.09222

    Returns:
    --------
    tuple: (cov_matrix_no_lumi, lumi_cov)
    cov_matrix_no_lumi: np.array, the sum of stat and syst covmats -> to be decomposed into artificial systematics
    lumi_cov: np.array, the lumi covmat. This is correlated between experiments so needs to be saved with type: SPECIAL
    """
    corr_matrix = get_correlation_matrix()
    uncertainties = get_uncertainties()

    # build correlated systematics covariance
    sys = np.array([uncertainties[i, 1] for i in range(3)])
    cov_sys = corr_matrix * np.outer(sys, sys)

    # array of lumi uncertainties
    lumi_unc = np.array([uncertainties[i, 2] for i in range(3)])

    # array of stat uncertainties
    stat = np.array([uncertainties[i, 0] for i in range(3)])

    return stat, cov_sys, lumi_unc


def get_systematics():
    stat, cov_sys, lumi_unc = get_covariance_matrices()

    # decompose sys covmat
    syst_unc = decompose_covmat(cov_sys)

    uncertainties = []

    # store only systematics for W+ and W-
    for i in range(3):
        uncertainties.append([{"name": f"ATLAS_WZ_TOT_13TEV_{i}", "values": [syst_unc[2, i]]}])

    uncertainties.append([{"name": "stat", "values": [stat[2]]}])
    uncertainties.append([{"name": "ATLAS_LUMI", "values": [lumi_unc[2]]}])

    return uncertainties
