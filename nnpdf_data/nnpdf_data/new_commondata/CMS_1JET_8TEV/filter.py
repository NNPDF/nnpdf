"""
Filter for CMS_1JET_8TEV

Created on Apr  2023

@author: Mark N. Costantini
"""

from filter_utils import (
    block_diagonal_corr,
    correlation_to_covariance,
    decompose_covmat,
    get_data_values,
    get_kinematics,
    get_stat_uncertainties,
    process_err,
    uncertainties_df,
)
import numpy as np
import yaml


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

    # construct block diagonal statistical covariance matrix (stat correlations within the same rapidity bins)
    stat_unc = get_stat_uncertainties()  # df_unc['ignore'].values
    # stat_unc = df_unc['stat'].values * df_unc['Sigma'].values / 100

    bd_stat_cov = correlation_to_covariance(block_diagonal_corr(tables), stat_unc)
    # bd_stat_cov = np.diag(stat_unc**2)

    # generate artificial systematics by decomposing statistical covariance matrix
    A_art_stat = decompose_covmat(bd_stat_cov)
    A_art_stat = np.nan_to_num(A_art_stat)  # set nan to zero

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

    # # NP corrections (SKIP)
    # np_p = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr+'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np = np.einsum('i,j->ij', np_p, np_p)
    # np_m = df_unc['Sigma'].values * (df_unc['NPCorr'].values * (1. + df_unc['npcorerr-'].values / 100.) - 1) / np.sqrt(2.)
    # cov_np += np.einsum('i,j->ij', np_m, np_m)

    covmat = cov_JES + cov_unfold + bd_stat_cov + lumi_cov + cov_uncorr

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

    return covmat


if __name__ == "__main__":
    # write data central values and kinematics to file
    filter_CMS_1JET_8TEV_data_kinetic()
    filter_CMS_1JET_8TEV_uncertainties()

    covmat = filter_CMS_1JET_8TEV_uncertainties()

    from validphys.api import API

    # why does the API give a 185 x 185 shaped covmat
    # i.e. why is it ignoring the low pt stuff

    inps = [{'dataset': "CMS_1JET_8TEV"}]
    inp = dict(dataset_inputs=inps, theoryid=200, use_cuts="internal")
    cmat = API.dataset_inputs_covmat_from_systematics(**inp)

    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Create a shared colorbar axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

    sns.heatmap(
        covmat / np.outer(np.sqrt(np.diag(covmat)), np.sqrt(np.diag(covmat))),
        annot=False,
        cmap="YlGnBu",
        ax=axs[0],
        cbar_ax=cbar_ax,
    )
    sns.heatmap(
        cmat / np.outer(np.sqrt(np.diag(cmat)), np.sqrt(np.diag(cmat))),
        annot=False,
        cmap="YlGnBu",
        ax=axs[1],
        cbar_ax=cbar_ax,
    )

    plt.show()

    ones = covmat / cmat
    print(ones)
    print()
    print(np.diag(ones))
    print()
    print(np.allclose(np.ones(covmat.shape), ones, rtol=1e-5))
    print(np.max(covmat - cmat), np.min(covmat - cmat))
    # print(np.argmax(covmat-cmat, axis=1), np.argmax(covmat-cmat, axis=0), np.argmin(covmat-cmat))
    max_index = np.argmax(covmat - cmat)
    print("Index of the largest entry:", max_index)

    # Get the row and column indices of the largest entry
    max_row_index, max_col_index = np.unravel_index(max_index, (covmat - cmat).shape)
    print("Row index of the largest entry:", max_row_index)
    print("Column index of the largest entry:", max_col_index)
    print((covmat.flatten()[max_index]), covmat[max_row_index, max_col_index])
    print((cmat.flatten()[max_index]), cmat[max_row_index, max_col_index])
    # print(covmat[np.argmax(covmat-cmat)])
