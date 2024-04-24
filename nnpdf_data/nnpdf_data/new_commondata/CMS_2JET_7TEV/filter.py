"""
Filter for CMS_2JET_7TEV

Created on Mar  2023
"""

import numpy as np
import pandas as pd
from scipy.linalg import block_diag
import yaml

from nnpdf_data.filter_utils.legacy_jets_utils import (
    JEC_error_matrix_CMS_2JET_7TEV,
    bin_by_bin_covmat_CMS_2JET_7TEV,
    dat_file_to_df_CMS_2JET_7TEV,
    get_corr_dat_file_CMS_2JET_7TEV,
    get_data_values_CMS_2JET_7TEV,
    get_kinematics_CMS_2JET_7TEV,
    get_stat_uncertainties_CMS_2JET_7TEV,
    lumi_covmat_CMS_2JET_7TEV,
    unfolding_error_matrix_CMS_2JET_7TEV,
)
from nnpdf_data.filter_utils.utils import correlation_to_covariance, decompose_covmat


def filter_CMS_2JET_7TEV_data_kinetic():
    """
    writes kinetic and data central values
    to kinematics.yaml and data.yaml files
    respectively
    """

    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get kinematics from hepdata tables
    kin = get_kinematics_CMS_2JET_7TEV(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values_CMS_2JET_7TEV(tables, version)

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filterCMS_2JET_7TEV_uncertainties():
    """
    writes uncertainties to uncertainties.yaml file
    """
    # read metadata file
    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['hepdata']['tables']
    version = metadata['hepdata']['version']

    # generate block diagonal statistical covariance matrix
    # Statistical uncertainty correlated between the mass bins in
    # the same rapidity range

    # get correlation matrix for statistical uncertainties
    corr_matrices = get_corr_dat_file_CMS_2JET_7TEV('rawdata/dijet_corr.dat')
    # get statistical uncertainties from each HEPData table
    stat_uncertainties = get_stat_uncertainties_CMS_2JET_7TEV()

    stat_cov_mats = []

    for i, table in enumerate(tables):
        if corr_matrices[i].shape[0] != np.array(stat_uncertainties[table]).shape[0]:
            raise ("Shapes of correlation matrix and uncertainties array are not compatible")
        # convert correlation matrices to covariance matrices
        stat_cov_mats.append(correlation_to_covariance(corr_matrices[i], stat_uncertainties[table]))

    # build block diagonal stat covmat
    BD_stat = stat_cov_mats[0]
    for i in range(1, len(stat_cov_mats)):
        stat = stat_cov_mats[i]
        BD_stat = block_diag(BD_stat, stat)

    # dataframe of uncertainties
    dfs = dat_file_to_df_CMS_2JET_7TEV()
    df_uncertainties = pd.concat(dfs, axis=0)
    cv = get_data_values_CMS_2JET_7TEV(tables, version)
    cv = np.array(cv)

    # get Luminosity Covmat CMSLUMI11
    lumi_cov = lumi_covmat_CMS_2JET_7TEV()
    A_lum = df_uncertainties["Lumi+"].multiply(cv, axis=0).to_numpy()

    # Get JEC covmat, CORR
    jec_error_matrix = JEC_error_matrix_CMS_2JET_7TEV()

    # get unfolding covmat, CORR
    unfold_error_matrix = unfolding_error_matrix_CMS_2JET_7TEV()

    # get bin-by-bin covmat, UNCORR
    bin_cov = bin_by_bin_covmat_CMS_2JET_7TEV()
    A_bin = df_uncertainties["Bin-by-bin-"].multiply(cv, axis=0).to_numpy()

    # generate artificial systematics
    A_art_sys_corr = decompose_covmat(covmat=BD_stat)

    # error definition
    error_definition = {
        f"art_sys_corr_{i}": {
            "description": f"artificial systematic originating from correlated statistical uncertainties",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(1, A_art_sys_corr.shape[0] + 1)
    }

    for col in jec_error_matrix.columns:
        error_definition[f"{col}"] = {
            "description": f"JEC uncertainty, {col}",
            "treatment": "MULT",
            "type": "CORR",
        }

    for col in unfold_error_matrix.columns:
        error_definition[f"{col}"] = {
            "description": f"Unfolding uncertainty, {col}",
            "treatment": "MULT",
            "type": "CORR",
        }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "MULT",
        "type": "CMSLUMI11",
    }

    error_definition["bin_by_bin_uncertainty"] = {
        "description": "bin_by_bin_uncertainty",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    # store error in dict
    error = []
    for n in range(A_art_sys_corr.shape[0]):
        error_value = {}

        for m in range(A_art_sys_corr.shape[1]):
            error_value[f"art_sys_corr_{m+1}"] = float(A_art_sys_corr[n, m])

        for col, m in zip(
            unfold_error_matrix.columns, range(unfold_error_matrix.to_numpy().shape[1])
        ):
            error_value[f"{col}"] = float(unfold_error_matrix.to_numpy()[n, m])

        for col, m in zip(jec_error_matrix.columns, range(jec_error_matrix.to_numpy().shape[1])):
            error_value[f"{col}"] = float(jec_error_matrix.to_numpy()[n, m])

        error_value["luminosity_uncertainty"] = float(A_lum[n])
        error_value["bin_by_bin_uncertainty"] = float(A_bin[n])

        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open(f"uncertainties.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # save data central values and kinematics
    filter_CMS_2JET_7TEV_data_kinetic()

    # save uncertainties
    filterCMS_2JET_7TEV_uncertainties()
