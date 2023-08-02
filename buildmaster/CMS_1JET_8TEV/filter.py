"""
Filter for CMS_1JET_8TEV

Created on Apr  2023

@author: Mark N. Costantini
"""

import yaml
import numpy as np

from filter_utils import (
    get_data_values,
    get_kinematics,
    get_stat_uncertainties,
    block_diagonal_corr,
    correlation_to_covariance,
    uncertainties_df,
    process_err,
)


def filter_CMS_1JET_8TEV_data_kinetic():
    """
    writes kinetic and data central values
    to kinematics.yaml and data.yaml files
    respectively
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get kinematics from hepdata tables
    kin = get_kinematics(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values(tables, version)

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_CMS_1JET_8TEV_uncertainties():
    """ """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get dataframe of uncertainties
    df_unc = uncertainties_df(tables)

    # construct block diagonal statistical covariance matrix
    stat_unc = df_unc['ignore'].values
    # stat_unc = df_unc['stat'].values * df_unc['Sigma'].values / 100

    bd_stat_cov = correlation_to_covariance(
        block_diagonal_corr(tables), stat_unc
    )  # get_stat_uncertainties())

    # Luminosity uncertainty
    lum_unc = df_unc['Luminosity'].values * df_unc['Sigma'].values / 100
    lumi_cov = np.outer(lum_unc, lum_unc)  # np.einsum('i,j->ij', lum_unc, lum_unc)

    # Uncorrelated
    uncorr = df_unc['uncor'].values * df_unc['Sigma'].values / 100
    cov_uncorr = np.diag(uncorr**2)

    # Unfolding Systematics
    df_unfold = (
        df_unc[['Unfolding+', 'Unfolding-']]
        * df_unc['Sigma'].values[:, np.newaxis]
        / 100.0
        / np.sqrt(2.0)
    )
    cov_unfold = df_unfold.to_numpy() @ df_unfold.to_numpy().T

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

    df_JES = process_err(df_JES)

    cov_JES = df_JES.to_numpy() @ df_JES.to_numpy().T

    import IPython

    IPython.embed()

    # # NP corrections (SKIP)
    # np_p = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr+'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np = np.einsum('i,j->ij', np_p, np_p)
    # np_m = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr-'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np += np.einsum('i,j->ij', np_m, np_m)

    covmat = cov_JES + cov_unfold + bd_stat_cov + lumi_cov + cov_uncorr

    return covmat


if __name__ == "__main__":
    # write data central values and kinematics to file
    # filter_CMS_1JET_8TEV_data_kinetic()
    covmat = filter_CMS_1JET_8TEV_uncertainties()

    from validphys.api import API

    # why does the API give a 185 x 185 shaped covmat
    # i.e. why is it ignoring the low pt stuff

    inps = [{'dataset': "CMS_1JET_8TEV"}]
    inp = dict(dataset_inputs=inps, theoryid=200, use_cuts="internal")
    cmat = API.dataset_inputs_covmat_from_systematics(**inp)

    print(cmat / covmat)
    print()
    print(np.diag(covmat / cmat))
    print()
    print(np.allclose(np.ones(covmat.shape), covmat / cmat))
