"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np
from sklearn.decomposition import PCA

from validphys import covmats
from validphys.calcutils import calc_chi2
from validphys.results import ThPredictionsResult

from reportengine import collect


def principal_components_dataset(dataset, fits_pdf, variancepdf, explained_variance_ratio=0.99):
    """
    Compute the principal components of theory predictions replica matrix
    (Ndat x Nrep feature matrix).

    Parameters
    ----------
    dataset: (DataSetSpec, DataGroupSpec)
        dataset for which the theory predictions and t0 covariance matrix
        will be loaded. Note that due to the structure of `validphys` this
        function can be overloaded to accept a DataGroupSpec.

    fits_pdf: list
        list of PDF objects produced from performing multiple closure tests
        fits. Each fit should have a different filterseed but the same
        underlying law used to generate the pseudodata.

    variancepdf: validphys.core.PDF
            PDF object used to estimate the variance of the fits.

    explained_variance_ratio: float, default is 0.93

    Returns
    -------
    tuple
        3D tuple:
        - matrix of the principal components (PCs) of shape (N_pc, N_dat)
        - reduced feature matrix, i.e., feature matrix projected onto PCs of shape (N_pc, N_rep)
        - N_pc: number of principal components kept

    """
    # get replicas from variance fit, used to estimate variance
    reps = ThPredictionsResult.from_convolution(variancepdf, dataset).error_members

    # something that could be tested: rescale feature matrix
    # reps_scaled = reps.preprocessing.scale(reps)

    # choose number of principal components (PCs) based on explained variance ratio
    n_comp = 1
    for _ in range(reps.shape[0]):
        pca = PCA(n_comp).fit(reps.T)
        if np.sum(pca.explained_variance_ratio_) >= explained_variance_ratio:
            break
        n_comp += 1

    # project feature matrix onto PCs
    pc_reps = pca.components_ @ reps

    return pca.components_, pc_reps, n_comp


def principal_components_bias_variance_dataset(
    internal_multiclosure_dataset_loader, explained_variance_ratio=0.99
):
    """
    Compute the bias and variance for one datasets
    using the principal component reduced covariance matrix.

    Parameters
    ----------
    internal_multiclosure_dataset_loader : tuple
        Tuple containing the results of multiclosure fits

    principal_components_dataset : tuple
        3D tuple containing the principal components of the theory predictions

    Returns
    -------
    tuple
        3D tuple:
        - biases: 1-D array of shape (Nfits,)
        - variances: 1-D array of shape (Nfits, )
        - n_comp: number of principal components kept
    """

    closures_th, law_th, _, _ = internal_multiclosure_dataset_loader

    reps = np.asarray([th.error_members for th in closures_th])

    # compute the covariance matrix of the theory predictions for each fit
    _covmats = [np.cov(rep, rowvar=True) for rep in reps]

    # compute the mean covariance matrix
    _covmat_mean = np.mean(_covmats, axis=0)

    # diagonalize the mean covariance matrix and only keep the principal components
    # that explain the required variance

    if _covmat_mean.shape == ():
        return np.nan, np.nan, 1

    eighvals, eigvecs = np.linalg.eigh(_covmat_mean)
    idx = np.argsort(eighvals)[::-1]
    # sort eigenvalues from largest to smallest
    eigvecs = eigvecs[:, idx]
    eighvals = eighvals[idx]
    eighvals_norm = eighvals / eighvals.sum()

    # choose components to keep based on EVR
    n_comp = 1
    for _ in range(eighvals.shape[0]):
        if np.sum(eighvals_norm[:n_comp]) >= explained_variance_ratio:
            break
        n_comp += 1

    # get the principal components
    pc_basis = eigvecs[:, :n_comp]

    # compute the (PCA) regularized covariance matrix
    covmat_pdf = pc_basis.T @ _covmat_mean @ pc_basis

    if n_comp <= 1:
        return np.nan, np.nan, n_comp

    # compute sqrt of pdf covariance matrix
    sqrt_covmat_pdf = covmats.sqrt_covmat(covmat_pdf)

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - law_th.central_value[:, np.newaxis]
    # shape here is (Npc, Nfits)
    delta_bias = pc_basis.T @ delta_bias

    # compute biases, shape of biases is (Nfits)
    biases = calc_chi2(sqrt_covmat_pdf, delta_bias)

    variances = []
    for i in range(reps.shape[0]):
        diffs = pc_basis.T @ (reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True))
        variances.append(np.mean(calc_chi2(sqrt_covmat_pdf, diffs)))

    return biases, np.asarray(variances), n_comp


principal_components_bias_variance_datasets = collect(
    "principal_components_bias_variance_dataset", ("data",)
)
