"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np
import dataclasses
import pandas as pd

from validphys import covmats
from validphys.calcutils import calc_chi2
from validphys.results import ThPredictionsResult
from validphys.closuretest.closure_checks import check_multifit_replicas
from validphys.closuretest.multiclosure import bootstrapped_internal_multiclosure_dataset_loader

from reportengine import collect


@dataclasses.dataclass(frozen=True)
class PCAInternalMulticlosureLoader:
    """
    Parameters
    ----------
    closures_th: list
        list containing validphys.results.ThPredictionsResult objects
        for each fit

    law_th: ThPredictionsResult object
        underlying law theory predictions

    pc_basis: np.array
        basis of principal components

    n_comp: int
        number of principal components kept after regularisation

    covmat_pca: np.array
        regularised covariance matrix computed from replicas
        of theory predictions

    sqrt_covmat_pca: np.array
        cholesky decomposed covariance matrix
    """

    closures_th: list
    law_th: ThPredictionsResult
    pc_basis: np.array
    n_comp: int
    covmat_pca: np.array
    sqrt_covmat_pca: np.array


@check_multifit_replicas
def internal_multiclosure_dataset_loader_pca(
    internal_multiclosure_dataset_loader,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Similar to multiclosure.internal_multiclosure_dataset_loader but returns
    PCA regularised covariance matrix, where the covariance matrix has been computed
    from the replicas of the theory predictions.

    Parameters
    ----------
    internal_multiclosure_dataset_loader: tuple
        closure fits theory predictions,
        underlying law theory predictions,
        covariance matrix,
        sqrt covariance matrix

    explained_variance_ratio: float, default is 0.99

    _internal_max_reps: int, default is None
        Maximum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    _internal_min_reps: int, default is 20
        Minimum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    Returns
    -------
    PCAInternalMulticlosureLoader
    """
    closures_th, law_th, _, _ = internal_multiclosure_dataset_loader

    reps = np.asarray([th.error_members for th in closures_th])

    # compute the covariance matrix of the theory predictions for each fit
    _covmats = [np.cov(rep, rowvar=True, bias=True) for rep in reps]

    # compute the mean covariance matrix
    _covmat_mean = np.mean(_covmats, axis=0)

    # diagonalize the mean covariance matrix and only keep the principal components
    # that explain the required variance

    if _covmat_mean.shape == ():
        return PCAInternalMulticlosureLoader(
            closures_th=closures_th,
            law_th=law_th,
            pc_basis=1,
            n_comp=1,
            covmat_pca=_covmat_mean,
            sqrt_covmat_pca=np.sqrt(_covmat_mean),
        )

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
    covmat_pca = pc_basis.T @ _covmat_mean @ pc_basis

    if n_comp == 1:
        return PCAInternalMulticlosureLoader(
            closures_th=closures_th,
            law_th=law_th,
            pc_basis=pc_basis,
            n_comp=1,
            covmat_pca=covmat_pca,
            sqrt_covmat_pca=np.sqrt(covmat_pca),
        )

    # compute sqrt of pdf covariance matrix
    sqrt_covmat_pca = covmats.sqrt_covmat(covmat_pca)

    return PCAInternalMulticlosureLoader(
        closures_th=closures_th,
        law_th=law_th,
        pc_basis=pc_basis,
        n_comp=n_comp,
        covmat_pca=covmat_pca,
        sqrt_covmat_pca=sqrt_covmat_pca,
    )
    
@check_multifit_replicas
def internal_multiclosure_data_loader_pca(
    internal_multiclosure_data_loader,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Like multiclosure.internal_multiclosure_dataset_loader_pca except for all data

    Parameters
    ----------
    internal_multiclosure_data_loader: tuple
        closure fits theory predictions,
        underlying law theory predictions,
        covariance matrix,
        sqrt covariance matrix

    explained_variance_ratio: float, default is 0.99

    _internal_max_reps: int, default is None
        Maximum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    _internal_min_reps: int, default is 20
        Minimum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    Returns
    -------
    PCAInternalMulticlosureLoader
    """
    closures_th, law_th, _, _ = internal_multiclosure_data_loader
    reps = np.asarray([th.error_members for th in closures_th])
    
    # compute the covariance matrix of the theory predictions for each fit
    _covmats = [np.cov(rep, rowvar=True, bias=True) for rep in reps]
    # compute the mean covariance matrix
    _covmat_mean = np.mean(_covmats, axis=0)
    # Keep the sqrt of the diagonals to reconstruct the covmat later
    D = np.sqrt(np.diag(_covmat_mean))
    
    # compute the correlation matrix of the theory predictions for each fit
    _corrmats = [np.corrcoef(rep, rowvar=True, bias=True) for rep in reps]
    _corrmat_mean = np.mean(_corrmats, axis=0)
    # diagonalize the mean correlation matrix and only keep the principal components
    # that explain the required variance

    if _covmat_mean.shape == ():
        return PCAInternalMulticlosureLoader(
            closures_th=closures_th,
            law_th=law_th,
            pc_basis=1,
            n_comp=1,
            covmat_pca=_covmat_mean,
            sqrt_covmat_pca=np.sqrt(_covmat_mean),
        )

    eighvals, eigvecs = np.linalg.eigh(_corrmat_mean)
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

    # compute the (PCA) regularized correlation matrix
    corrmat_pca = pc_basis.T @ _corrmat_mean @ pc_basis
    # project the diagonal matrix into the PCA space 
    D_pca = pc_basis.T @ np.diag(D) @ pc_basis
    # compute the (PCA) regularized covariance matrix
    covmat_pca = D_pca @ corrmat_pca @ D_pca
    if n_comp == 1:
        return PCAInternalMulticlosureLoader(
            closures_th=closures_th,
            law_th=law_th,
            pc_basis=pc_basis,
            n_comp=1,
            covmat_pca=covmat_pca,
            sqrt_covmat_pca=np.sqrt(covmat_pca),
        )

    # compute sqrt of pdf covariance matrix
    sqrt_covmat_pca = covmats.sqrt_covmat(covmat_pca)
    return PCAInternalMulticlosureLoader(
        closures_th=closures_th,
        law_th=law_th,
        pc_basis=pc_basis,
        n_comp=n_comp,
        covmat_pca=covmat_pca,
        sqrt_covmat_pca=sqrt_covmat_pca,
    )


def bootstrapped_internal_multiclosure_dataset_loader_pca(
    internal_multiclosure_dataset_loader,
    n_fit_max,
    n_fit,
    n_rep_max,
    n_rep,
    n_boot_multiclosure,
    rng_seed_mct_boot,
    use_repeats=True,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Similar to multiclosure.bootstrapped_internal_multiclosure_dataset_loader but returns
    PCA regularised covariance matrix, where the covariance matrix has been computed
    from the replicas of the theory predictions.
    """

    # get bootstrapped internal multiclosure dataset loader
    bootstrap_imdl = bootstrapped_internal_multiclosure_dataset_loader(
        internal_multiclosure_dataset_loader,
        n_fit_max=n_fit_max,
        n_fit=n_fit,
        n_rep_max=n_rep_max,
        n_rep=n_rep,
        n_boot_multiclosure=n_boot_multiclosure,
        rng_seed_mct_boot=rng_seed_mct_boot,
        use_repeats=use_repeats,
    )

    # PCA regularise all the bootstrapped internal multiclosure dataset loaders
    bootstrap_imdl_pca = [
        internal_multiclosure_dataset_loader_pca(
            imdl, explained_variance_ratio, _internal_max_reps, _internal_min_reps
        )
        for imdl in bootstrap_imdl
    ]
    return tuple(bootstrap_imdl_pca)


def bootstrapped_internal_multiclosure_data_loader_pca(
    internal_multiclosure_data_loader,
    n_fit_max,
    n_fit,
    n_rep_max,
    n_rep,
    n_boot_multiclosure,
    rng_seed_mct_boot,
    use_repeats=True,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Same as bootstrapped_internal_multiclosure_dataset_loader_pca but for all the data.
    """
    return bootstrapped_internal_multiclosure_dataset_loader_pca(
        internal_multiclosure_data_loader,
        n_fit_max,
        n_fit,
        n_rep_max,
        n_rep,
        n_boot_multiclosure,
        rng_seed_mct_boot,
        use_repeats,
        explained_variance_ratio,
        _internal_max_reps,
        _internal_min_reps,
    )


def principal_components_bias_variance_dataset(internal_multiclosure_dataset_loader_pca):
    """
    Compute the bias and variance for one datasets
    using the principal component reduced covariance matrix.

    Parameters
    ----------
    internal_multiclosure_dataset_loader : tuple
        Tuple containing the results of multiclosure fits

    explained_variance_ratio : float, default is 0.99
        3D tuple containing the principal components of the theory predictions

    Returns
    -------
    tuple
        3D tuple:
        - biases: 1-D array of shape (Nfits,)
        - variances: 1-D array of shape (Nfits, )
        - n_comp: number of principal components kept
    """

    pca_loader = internal_multiclosure_dataset_loader_pca

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]

    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        biases = (delta_bias / pca_loader.sqrt_covmat_pca) ** 2
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean((diffs / pca_loader.sqrt_covmat_pca) ** 2))
    else:
        delta_bias = pca_loader.pc_basis.T @ delta_bias
        biases = calc_chi2(pca_loader.sqrt_covmat_pca, delta_bias)
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis.T @ (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean(calc_chi2(pca_loader.sqrt_covmat_pca, diffs)))

    return biases, np.asarray(variances), pca_loader.n_comp

def principal_components_bias_variance_data(internal_multiclosure_data_loader_pca):
    """
    Like principal_components_bias_variance_datasets but for all data

    Parameters
    ----------
    internal_multiclosure_data_loader_pca : tuple
        Tuple containing the results of multiclosure fits after pca regularization


    Returns
    -------
    tuple
        3D tuple:
        - biases: 1-D array of shape (Nfits,)
        - variances: 1-D array of shape (Nfits, )
        - n_comp: number of principal components kept
    """

    pca_loader = internal_multiclosure_data_loader_pca

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]

    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        biases = (delta_bias / pca_loader.sqrt_covmat_pca) ** 2
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean((diffs / pca_loader.sqrt_covmat_pca) ** 2))
    else:
        delta_bias = pca_loader.pc_basis.T @ delta_bias
        biases = calc_chi2(pca_loader.sqrt_covmat_pca, delta_bias)
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis.T @ (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean(calc_chi2(pca_loader.sqrt_covmat_pca, diffs)))

    return biases, np.asarray(variances), pca_loader.n_comp

def principal_components_normalized_delta_data(internal_multiclosure_data_loader_pca):
    """
    Compute for all data only the normalized delta after PCA regularization

    Parameters
    ----------
    internal_multiclosure_data_loader_pca : tuple
        Tuple containing the results of multiclosure fits after pca regularization


    Returns
    -------
    nd.array: deltas
    """

    pca_loader = internal_multiclosure_data_loader_pca

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]
    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        standard_deviations = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            standard_deviations.append(np.sqrt(np.mean((diffs / pca_loader.sqrt_covmat_pca) ** 2)))
    else:
        delta_bias = pca_loader.pc_basis.T @ delta_bias
        standard_deviations = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis.T @ (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            standard_deviations.append(np.std(diffs,axis=1))
    

    return (delta_bias/np.asarray(standard_deviations).T).flatten(), pca_loader.n_comp


principal_components_bias_variance_datasets = collect(
    "principal_components_bias_variance_dataset", ("data",)
)


def bootstrapped_principal_components_bias_variance_dataset(
    bootstrapped_internal_multiclosure_dataset_loader_pca, dataset
):
    """
    Computes Bias and Variance for each bootstrap sample.
    Returns a DataFrame with the results.
    """
    boot_bias_var_samples = []
    for i, boot_imdl_pca in enumerate(bootstrapped_internal_multiclosure_dataset_loader_pca):
        bias, var, n_comp = principal_components_bias_variance_dataset(boot_imdl_pca)
        boot_bias_var_samples.append(
            {
                "bias": np.mean(bias),
                "variance": np.mean(var),
                "n_comp": n_comp,
                "dataset": str(dataset),
                "bootstrap_index": i,
            }
        )

    df = pd.DataFrame.from_records(
        boot_bias_var_samples,
        index="bootstrap_index",
        columns=("bootstrap_index", "dataset", "n_comp", "bias", "variance"),
    )

    df.columns = ["dataset", "n_comp", "bias", "variance"]
    return df


bootstrapped_principal_components_bias_variance_data = collect(
    "bootstrapped_principal_components_bias_variance_dataset", ("data",)
)
