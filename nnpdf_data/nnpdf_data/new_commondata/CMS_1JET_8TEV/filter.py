"""
Filter for CMS_1JET_8TEV

Created on Apr  2023
"""

import numpy as np
import yaml

from nnpdf_data.filter_utils.legacy_jets_utils import (
    block_diagonal_corr_CMS_1JET_8TEV,
    get_data_values_CMS_1JET_8TEV,
    get_kinematics_CMS_1JET_8TEV,
    get_stat_uncertainties_CMS_1JET_8TEV,
    process_err_CMS_1JET_8TEV,
    uncertainties_df_CMS_1JET_8TEV,
)
from nnpdf_data.filter_utils.utils import correlation_to_covariance, decompose_covmat, prettify_float

yaml.add_representer(float, prettify_float)


def filter_CMS_1JET_8TEV_data_kinetic():
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
    kin = get_kinematics_CMS_1JET_8TEV(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values_CMS_1JET_8TEV(tables, version)

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_CMS_1JET_8TEV_uncertainties():
    """ """

    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get dataframe of uncertainties
    df_unc = uncertainties_df_CMS_1JET_8TEV(tables)

    # construct block diagonal statistical covariance matrix (stat correlations within the same rapidity bins)
    stat_unc = get_stat_uncertainties_CMS_1JET_8TEV()  # df_unc['ignore'].values
    # stat_unc = df_unc['stat'].values * df_unc['Sigma'].values / 100

    bd_stat_cov = correlation_to_covariance(block_diagonal_corr_CMS_1JET_8TEV(tables), stat_unc)
    # bd_stat_cov = np.diag(stat_unc**2)

    # generate artificial systematics by decomposing statistical covariance matrix
    A_art_stat = decompose_covmat(bd_stat_cov)
    A_art_stat = np.nan_to_num(A_art_stat)  # set nan to zero

    # Luminosity uncertainty
    lum_unc = df_unc['Luminosity'].values * df_unc['Sigma'].values / 100

    # Uncorrelated
    uncorr = df_unc['uncor'].values * df_unc['Sigma'].values / 100

    # Unfolding Systematics
    df_unfold = (
        df_unc[['Unfolding+', 'Unfolding-']]
        * df_unc['Sigma'].values[:, np.newaxis]
        / 100.0
        / np.sqrt(2.0)
    )

    # Systematic Correlated Uncertainties (Unfolding + JES)
    df_JES = (
        df_unc.iloc[:, 10:].drop(
            ['stat', 'uncor', 'Luminosity', 'Unfolding+', 'Unfolding-'], axis=1
        )
        * df_unc['Sigma'].values[:, np.newaxis]
        / 100.0
        / np.sqrt(2.0)
    )
    # if a pair of uncertainties has the same sign keep only the one with the
    # largest absolute value and set the other one to zero.

    df_JES = process_err_CMS_1JET_8TEV(df_JES)

    # # NP corrections (SKIP)
    # np_p = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr+'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np = np.einsum('i,j->ij', np_p, np_p)
    # np_m = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr-'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np += np.einsum('i,j->ij', np_m, np_m)

    # save systematics to yaml file

    # error definition
    error_definition = {
        f"art_sys_{i}": {
            "description": f"artificial systematic {i}, generated from block diagonal statistical covariance matrix",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(1, A_art_stat.shape[0] + 1)
    }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "MULT",
        "type": "CMSLUMI12",
    }

    error_definition["uncorrelated_uncertainty"] = {
        "description": "uncorrelated systematic uncertainty",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    for col in df_unfold.columns:
        error_definition[f"{col}"] = {
            "description": f"correlated unfolding uncertainty, {col}",
            "treatment": "MULT",
            "type": "CORR",
        }

    for col in df_JES:
        error_definition[f"{col}"] = {
            "description": f"correlated JES uncertainty, {col}",
            "treatment": "MULT",
            "type": "CORR",
        }

    # store error in dict
    error = []
    for n in range(A_art_stat.shape[0]):
        error_value = {}

        # artificial stat uncertainties
        for m in range(A_art_stat.shape[1]):
            error_value[f"art_sys_{m+1}"] = float(A_art_stat[n, m])

        # unfolding uncertainties
        for col, m in zip(df_unfold.columns, range(df_unfold.to_numpy().shape[1])):
            error_value[f"{col}"] = float(df_unfold.to_numpy()[n, m])

        # JES uncertainties
        for col, m in zip(df_JES.columns, range(df_JES.to_numpy().shape[1])):
            error_value[f"{col}"] = float(df_JES.to_numpy()[n, m])

        # luminosity uncertainties
        error_value["luminosity_uncertainty"] = float(lum_unc[n])

        # uncorrelated uncertainties
        error_value["uncorrelated_uncertainty"] = float(uncorr[n])

        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open(f"uncertainties.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # write data central values and kinematics to file
    filter_CMS_1JET_8TEV_data_kinetic()
    filter_CMS_1JET_8TEV_uncertainties()
