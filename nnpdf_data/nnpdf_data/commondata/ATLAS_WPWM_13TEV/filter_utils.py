"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml
import numpy as np
from nnpdf_data.filter_utils.utils import decompose_covmat


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    mw2 = 80.385**2

    for i in range(2):

        kin_value = {
            'k1': {'min': None, 'mid': 0.0, 'max': None},
            'M2': {'min': None, 'mid': mw2, 'max': None},
            'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values():
    """
    returns the central data values in the form of a list.
    """
    hepdata_table_wp, hepdata_table_wm = (
        "rawdata/HEPData-ins1436497-v1-Table_8.yaml",
        "rawdata/HEPData-ins1436497-v1-Table_9.yaml",
    )

    with open(hepdata_table_wp, 'r') as file:
        input_wp = yaml.safe_load(file)

    with open(hepdata_table_wm, 'r') as file:
        input_wm = yaml.safe_load(file)

    values_wm = input_wm['dependent_variables'][0]['values']
    values_wp = input_wp['dependent_variables'][0]['values']

    data_central = [values_wm[0]['value'] * 1000000, values_wp[0]['value'] * 1000000]

    return data_central


def get_correlation_matrix():
    """
    See extra material page: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/STDM-2015-03/tabaux_03.pdf

    Note that this does not include the normalisation uncertainty due to the luminosity.
    """

    correlation_matrix = np.ones((2, 2))
    correlation_matrix[0, 1] = 0.93
    correlation_matrix[1, 0] = correlation_matrix[0, 1]

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
    cv = get_data_values()

    stat_wm = 0.01 * cv[0]
    stat_wp = 0.01 * cv[1]

    syst_wm = 0.07 * cv[0]
    syst_wp = 0.09 * cv[1]

    lumi_wm = 0.10 * cv[0]
    lumi_wp = 0.07 * cv[1]

    stat_cov = np.diag([stat_wm**2, stat_wp**2])
    # lumi_cov = np.einsum("i,j->ij", np.array([lumi_wm, lumi_wp]), np.array([lumi_wm, lumi_wp]))
    lumi_unc = np.array([lumi_wm, lumi_wp])
    syst_cov = corr_matrix * np.outer(np.array([syst_wm, syst_wp]), np.array([syst_wm, syst_wp]))

    cov_matrix_no_lumi = stat_cov + syst_cov

    return cov_matrix_no_lumi, lumi_unc


def get_systematics():
    """
    Does cholesky decomposition of syst + stat covmat and returns uncertainties
    list with artificial sys + lumi uncertainties.
    """
    cov_matrix_no_lumi, lumi_unc = get_covariance_matrices()

    # decompose covmat
    syst_unc = decompose_covmat(cov_matrix_no_lumi)

    uncertainties = []

    uncertainties.append([{"name": "stat", "values": [syst_unc[0, 0], syst_unc[1, 0]]}])
    uncertainties.append([{"name": "sys1", "values": [syst_unc[0, 1], syst_unc[1, 1]]}])
    uncertainties.append([{"name": "ATLAS_LUMI", "values": [lumi_unc[0], lumi_unc[1]]}])

    return uncertainties
