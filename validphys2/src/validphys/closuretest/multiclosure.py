"""
closuretest/multiclosure.py

Module containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit. The actions
in this module are used to produce results which are plotted in
``multiclosure_output.py``

"""

import numpy as np
import dataclasses
import pandas as pd

from reportengine import collect
import validphys
from validphys.calcutils import calc_chi2
from validphys import covmats
from validphys.checks import check_use_t0
from validphys.closuretest.closure_checks import (
    check_at_least_10_fits,
    check_fits_areclosures,
    check_fits_different_filterseed,
    check_fits_underlying_law_match,
    check_multifit_replicas,
    check_t0pdfset_matches_multiclosure_law,
)
from validphys.results import ThPredictionsResult

# bootstrap seed default
DEFAULT_SEED = 9689372
# stepsize in fits/replicas to use for finite size bootstraps
SAMPLING_INTERVAL = 5


@dataclasses.dataclass(frozen=True)
class MulticlosureLoader:
    """
    Stores the basic information for a multiclosure study.

    Attributes
    ----------
    closure_theories: list
        List of validphys.results.ThPredictionsResult objects for each fit.

    law_theory: validphys.results.ThPredictionsResult
        ThPredictionsResult object for the underlying law.

    covmat_reps_mean: np.array
        Covariance matrix of the theory predictions averaged over fits.
    """

    closure_theories: list
    law_theory: ThPredictionsResult
    covmat_reps_mean: np.array


@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
@check_use_t0
@check_t0pdfset_matches_multiclosure_law
def multiclosure_dataset_loader(
    dataset: validphys.core.DataSetSpec,
    fits_pdf: list,
    multiclosure_underlyinglaw: validphys.core.PDF,
    t0set: validphys.core.PDF,
) -> MulticlosureLoader:
    """
    Internal function for loading multiple theory predictions and underlying law
    for a given dataset. This function is used to avoid memory issues when
    caching the load function of a group of datasets.

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
    multiclosure_underlyinglaw: PDF
        PDF used to generate the pseudodata which the closure tests fitted. This
        is inferred from the fit runcards.
    t0set: validphys.core.PDF
        t0 pdfset, is only used to check that the underlying law matches the t0set.


    Returns
    -------
    MulticlosureLoader
        A dataclass storing the theory predictions for the fits and the underlying law.

    Notes
    -----
    This function replicates behaviour found elsewhere in validphys, the reason
    for this is that due to the default caching behaviour one can run into
    memory issues when loading the theory predictions for the amount of fits
    typically used in these studies.
    """
    closure_theories = [ThPredictionsResult.from_convolution(pdf, dataset) for pdf in fits_pdf]
    law_theory = ThPredictionsResult.from_convolution(multiclosure_underlyinglaw, dataset)

    reps = np.asarray([th.error_members for th in closure_theories])
    # compute the covariance matrix of the theory predictions for each fit
    _covmats = np.array([np.cov(rep, rowvar=True, bias=True) for rep in reps])

    # compute the mean covariance matrix
    covmat_reps_mean = np.mean(_covmats, axis=0)

    return MulticlosureLoader(
        closure_theories=closure_theories, law_theory=law_theory, covmat_reps_mean=covmat_reps_mean
    )


@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
@check_t0pdfset_matches_multiclosure_law
@check_use_t0
def multiclosure_data_loader(
    data: validphys.core.DataGroupSpec,
    fits_pdf: list,
    multiclosure_underlyinglaw: validphys.core.PDF,
    t0set: validphys.core.PDF,
) -> MulticlosureLoader:
    """Like `multiclosure_dataset_loader` except for all data"""
    return multiclosure_dataset_loader(data, fits_pdf, multiclosure_underlyinglaw, t0set)


def eigendecomposition(covmat: np.array) -> tuple:
    """
    Computes the eigendecomposition of a covariance matrix
    and returns the eigenvalues, eigenvectors and the normalized
    eigenvalues ordered from largest to smallest.

    Parameters
    ----------
    covmat: np.array
        covariance matrix

    Returns
    -------
    tuple
        3D tuple containing the eigenvalues, eigenvectors and the normalized
        eigenvalues.
        Note that the eigenvalues are sorted from largest to smallest.
    """
    eighvals, eigvecs = np.linalg.eigh(covmat)
    idx = np.argsort(eighvals)[::-1]
    # sort eigenvalues from largest to smallest
    eigvecs = eigvecs[:, idx]
    eighvals = eighvals[idx]
    eighvals_norm = eighvals / eighvals.sum()

    return eighvals, eigvecs, eighvals_norm


@dataclasses.dataclass(frozen=True)
class RegularizedMulticlosureLoader(MulticlosureLoader):
    """
    Attributes
    ----------
    pc_basis: np.array
        Basis of principal components.

    n_comp: int
        Number of principal components kept after regularisation.

    reg_covmat_reps_mean: np.array
        Diagonal, regularised covariance matrix computed from replicas
        of theory predictions.

    sqrt_reg_covmat_reps_mean: np.array
        Sqrt of the regularised covariance matrix.
    """

    pc_basis: np.array
    n_comp: int
    reg_covmat_reps_mean: np.array
    sqrt_reg_covmat_reps_mean: np.array


@check_multifit_replicas
def regularized_multiclosure_dataset_loader(
    multiclosure_dataset_loader: MulticlosureLoader,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
) -> RegularizedMulticlosureLoader:
    """
    Similar to multiclosure.multiclosure_dataset_loader but computes the regularized
    PDF covariance matrix by only keeping the largest eigenvalues that sum to the
    `explained_variance_ratio`.


    Parameters
    ----------
    multiclosure_dataset_loader: MulticlosureLoader

    explained_variance_ratio: float, default is 0.99

    _internal_max_reps: int, default is None
        Maximum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    _internal_min_reps: int, default is 20
        Minimum number of replicas used in the fits
        this is needed to check that the number of replicas is the same for all fits

    Returns
    -------
    RegularizedMulticlosureLoader
    """
    closures_th = multiclosure_dataset_loader.closure_theories
    law_th = multiclosure_dataset_loader.law_theory
    covmat_reps_mean = multiclosure_dataset_loader.covmat_reps_mean

    if covmat_reps_mean.shape == ():
        return RegularizedMulticlosureLoader(
            closure_theories=closures_th,
            law_theory=law_th,
            covmat_reps_mean=covmat_reps_mean,
            pc_basis=1,
            n_comp=1,
            reg_covmat_reps_mean=covmat_reps_mean,
            sqrt_reg_covmat_reps_mean=np.sqrt(covmat_reps_mean),
        )

    # diagonalize the mean covariance matrix and only keep the principal components
    # that explain the required variance
    eighvals, eigvecs, eighvals_norm = eigendecomposition(covmat_reps_mean)

    # choose components to keep based on EVR
    n_comp = 1
    for _ in range(eighvals.shape[0]):
        if np.sum(eighvals_norm[:n_comp]) >= explained_variance_ratio:
            break
        n_comp += 1

    # get the principal components
    pc_basis = eigvecs[:, :n_comp]

    # Diagonalise and project the mean covmat in the space spanned by the PCs
    reg_covmat_reps_mean = pc_basis.T @ covmat_reps_mean @ pc_basis

    if n_comp == 1:
        return RegularizedMulticlosureLoader(
            closure_theories=closures_th,
            law_theory=law_th,
            covmat_reps_mean=covmat_reps_mean,
            pc_basis=pc_basis,
            n_comp=1,
            reg_covmat_reps_mean=reg_covmat_reps_mean,
            sqrt_reg_covmat_reps_mean=np.sqrt(reg_covmat_reps_mean),
        )

    # compute sqrt of pdf covariance matrix (NOTE: should be the same as taking np.sqrt() since the matrix is diagonal)
    sqrt_reg_covmat_reps_mean = covmats.sqrt_covmat(reg_covmat_reps_mean)

    return RegularizedMulticlosureLoader(
        closure_theories=closures_th,
        law_theory=law_th,
        covmat_reps_mean=covmat_reps_mean,
        pc_basis=pc_basis,
        n_comp=n_comp,
        reg_covmat_reps_mean=reg_covmat_reps_mean,
        sqrt_reg_covmat_reps_mean=sqrt_reg_covmat_reps_mean,
    )


@check_multifit_replicas
def regularized_multiclosure_data_loader(
    multiclosure_data_loader,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Like multiclosure.regularized_multiclosure_dataset_loader except for all data

    Parameters
    ----------
    multiclosure_data_loader: tuple
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
    RegularizedMulticlosureLoader
    """
    closures_th = multiclosure_data_loader.closure_theories
    law_th = multiclosure_data_loader.law_theory
    covmat_reps_mean = multiclosure_data_loader.covmat_reps_mean

    # Keep the sqrt of the diagonals to reconstruct the covmat later
    D = np.sqrt(np.diag(covmat_reps_mean))

    # compute the correlation matrix
    _corrmat_mean = covmat_reps_mean / np.outer(D, D)

    # diagonalize the mean correlation matrix and only keep the principal components
    # that explain the required variance
    if covmat_reps_mean.shape == ():
        return RegularizedMulticlosureLoader(
            closure_theories=closures_th,
            law_theory=law_th,
            covmat_reps_mean=covmat_reps_mean,
            pc_basis=1,
            n_comp=1,
            reg_covmat_reps_mean=covmat_reps_mean,
            sqrt_reg_covmat_reps_mean=np.sqrt(covmat_reps_mean),
        )

    eighvals, eigvecs, eighvals_norm = eigendecomposition(_corrmat_mean)

    # choose components to keep based on EVR
    n_comp = 1
    for _ in range(eighvals.shape[0]):
        if np.sum(eighvals_norm[:n_comp]) >= explained_variance_ratio:
            break
        n_comp += 1
    # get the principal components
    pc_basis = eigvecs[:, :n_comp]

    # compute the (PCA) regularized covariance matrix by projecting the mean covariance matrix
    # onto the principal components
    reg_covmat_reps_mean = pc_basis.T @ covmat_reps_mean @ pc_basis

    if n_comp == 1:
        return RegularizedMulticlosureLoader(
            closure_theories=closures_th,
            law_theory=law_th,
            covmat_reps_mean=covmat_reps_mean,
            pc_basis=pc_basis,
            n_comp=1,
            reg_covmat_reps_mean=reg_covmat_reps_mean,
            sqrt_reg_covmat_reps_mean=np.sqrt(reg_covmat_reps_mean),
        )

    # compute sqrt of pdf covariance matrix
    sqrt_reg_covmat_reps_mean = covmats.sqrt_covmat(reg_covmat_reps_mean)
    return RegularizedMulticlosureLoader(
        closure_theories=closures_th,
        law_theory=law_th,
        covmat_reps_mean=covmat_reps_mean,
        pc_basis=pc_basis,
        n_comp=n_comp,
        reg_covmat_reps_mean=reg_covmat_reps_mean,
        sqrt_reg_covmat_reps_mean=sqrt_reg_covmat_reps_mean,
    )


def bootstrapped_internal_multiclosure_dataset_loader_pca(
    multiclosure_dataset_loader,
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
        multiclosure_dataset_loader,
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
        regularized_multiclosure_dataset_loader(
            imdl, explained_variance_ratio, _internal_max_reps, _internal_min_reps
        )
        for imdl in bootstrap_imdl
    ]
    return tuple(bootstrap_imdl_pca)


def bootstrapped_internal_multiclosure_data_loader_pca(
    multiclosure_data_loader,
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
    Same as `bootstrapped_internal_multiclosure_dataset_loader_pca` but for all the data.
    """
    # get bootstrapped internal multiclosure dataset loader
    bootstrap_imdl = bootstrapped_internal_multiclosure_data_loader(
        multiclosure_data_loader,
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
        regularized_multiclosure_data_loader(
            imdl, explained_variance_ratio, _internal_max_reps, _internal_min_reps
        )
        for imdl in bootstrap_imdl
    ]
    return tuple(bootstrap_imdl_pca)


def principal_components_bias_variance_dataset(regularized_multiclosure_dataset_loader):
    """
    Compute the bias and variance for one dataset
    using the principal component reduced covariance matrix.

    Parameters
    ----------
    multiclosure_dataset_loader : tuple
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

    pca_loader = regularized_multiclosure_dataset_loader

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]

    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        biases = (delta_bias / pca_loader.sqrt_reg_covmat_reps_mean) ** 2
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean((diffs / pca_loader.sqrt_reg_covmat_reps_mean) ** 2))
    else:
        delta_bias = pca_loader.pc_basis.T @ delta_bias
        biases = calc_chi2(pca_loader.sqrt_reg_covmat_reps_mean, delta_bias)
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis.T @ (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean(calc_chi2(pca_loader.sqrt_reg_covmat_reps_mean, diffs)))

    return biases, np.asarray(variances), pca_loader.n_comp


def principal_components_bias_variance_data(regularized_multiclosure_data_loader):
    """
    Like principal_components_bias_variance_datasets but for all data

    Parameters
    ----------
    regularized_multiclosure_data_loader : tuple
        Tuple containing the results of multiclosure fits after pca regularization


    Returns
    -------
    tuple
        3D tuple:
        - biases: 1-D array of shape (Nfits,)
        - variances: 1-D array of shape (Nfits, )
        - n_comp: number of principal components kept
    """

    pca_loader = regularized_multiclosure_data_loader

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]

    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        biases = (delta_bias / pca_loader.sqrt_reg_covmat_reps_mean) ** 2
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean((diffs / pca_loader.sqrt_reg_covmat_reps_mean) ** 2))
    else:
        delta_bias = pca_loader.pc_basis.T @ delta_bias
        biases = calc_chi2(pca_loader.sqrt_reg_covmat_reps_mean, delta_bias)
        variances = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis.T @ (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            variances.append(np.mean(calc_chi2(pca_loader.sqrt_reg_covmat_reps_mean, diffs)))

    return biases, np.asarray(variances), pca_loader.n_comp


def principal_components_normalized_delta_data(regularized_multiclosure_data_loader):
    """
    Compute for all data only the normalized delta after PCA regularization

    Parameters
    ----------
    regularized_multiclosure_data_loader : tuple
        Tuple containing the results of multiclosure fits after pca regularization


    Returns
    -------
    nd.array: deltas
    """

    pca_loader = regularized_multiclosure_data_loader

    reps = np.asarray([th.error_members for th in pca_loader.closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - pca_loader.law_th.central_value[:, np.newaxis]

    # find basis that diagonalise covmat pca
    eigvals, eigenvects = np.linalg.eigh(pca_loader.reg_covmat_reps_mean)

    if pca_loader.n_comp == 1:
        delta_bias = pca_loader.pc_basis * delta_bias
        std_deviations = []
        for i in range(reps.shape[0]):
            diffs = pca_loader.pc_basis * (
                reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
            )
            std_deviations.append(
                np.sqrt(np.mean((diffs / pca_loader.sqrt_reg_covmat_reps_mean) ** 2))
            )
    else:
        delta_bias = eigenvects.T @ (pca_loader.pc_basis.T @ delta_bias)
        std_deviations = np.sqrt(eigvals)[:, None]

    return (delta_bias / std_deviations).flatten(), pca_loader.n_comp


principal_components_bias_variance_datasets = collect(
    "principal_components_bias_variance_dataset", ("data",)
)


def bootstrapped_principal_components_normalized_delta_data(
    bootstrapped_internal_multiclosure_data_loader_pca,
):
    """
    Compute the normalized deltas for each bootstrap sample.

    Parameters
    ----------
    bootstrapped_internal_multiclosure_data_loader_pca: list
        list of tuples containing the results of multiclosure fits after pca regularization

    Returns
    -------
    list
        list of tuples containing the normalized deltas and the number of principal components.
        Each tuple corresponds to a bootstrap sample.
    """
    normalised_deltas = []
    for boot_imdl_pca in bootstrapped_internal_multiclosure_data_loader_pca:
        normalised_deltas.append(principal_components_normalized_delta_data(boot_imdl_pca))
    return normalised_deltas


def bootstrapped_indicator_function_data(
    bootstrapped_principal_components_normalized_delta_data, nsigma=1
):
    """
    Compute the indicator function for each bootstrap sample.

    Parameters
    ----------
    bootstrapped_principal_components_normalized_delta_data: list
        list of tuples containing the normalized deltas and the number of principal components.
        Each tuple corresponds to a bootstrap sample.

    nsigma: int, default is 1

    Returns
    -------
    2-D tuple:
        list
            list of length N_boot and entrances are arrays of dim Npca x Nfits containing the indicator function for each bootstrap sample.

        float
            average number of degrees of freedom
    """
    indicator_list = []
    ndof_list = []
    for boot, ndof in bootstrapped_principal_components_normalized_delta_data:
        indicator_list.append(standard_indicator_function(boot, nsigma))
        ndof_list.append(ndof)
    return indicator_list, np.mean(np.asarray(ndof_list))


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


bootstrapped_principal_components_bias_variance_datasets = collect(
    "bootstrapped_principal_components_bias_variance_dataset", ("data",)
)


def bootstrapped_principal_components_bias_variance_data(
    bootstrapped_internal_multiclosure_data_loader_pca,
):
    """
    Computes Bias and Variance for each bootstrap sample.
    Returns a DataFrame with the results.
    """
    boot_bias_var_samples = []
    for i, boot_imdl_pca in enumerate(bootstrapped_internal_multiclosure_data_loader_pca):
        bias, var, n_comp = principal_components_bias_variance_data(boot_imdl_pca)
        boot_bias_var_samples.append(
            {
                "bias": np.mean(bias),
                "variance": np.mean(var),
                "n_comp": n_comp,
                "data": "Full dataset",
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


@check_multifit_replicas
def fits_normed_dataset_central_delta(
    multiclosure_dataset_loader, _internal_max_reps=None, _internal_min_reps=20
):
    """
    For each fit calculate the difference between central expectation value and true val. Normalize this
    value by the variance of the differences between replicas and central expectation value (different
    for each fit but expected to vary only a little). Each observable central exp value is
    expected to be gaussianly distributed around the true value set by the fakepdf.

    Parameters
    ----------
    multiclosure_dataset_loader: tuple
        closure fits theory predictions,
        underlying law theory predictions,
        covariance matrix,
        sqrt covariance matrix

    _internal_max_reps: int
        maximum number of replicas to use for each fit

    _internal_min_reps: int
        minimum number of replicas to use for each fit

    Returns
    -------
    deltas: np.array
        2-D array with shape (n_fits, n_obs)
    """
    closures_th = multiclosure_dataset_loader.closure_theories
    law_th = multiclosure_dataset_loader.law_theory

    # The dimentions here are (fit, data point, replica)
    reps = np.asarray([th.error_members[:, :_internal_max_reps] for th in closures_th])
    # One could mask here some reps in order to avoid redundancy of information
    # TODO

    n_fits = np.shape(reps)[0]
    deltas = []
    # There are n_fits pdf_covariances
    # flag to see whether to eliminate dataset
    for rep in reps:
        # bias diffs in the for loop should have shape (n_obs,)
        bias_diffs = np.mean(rep, axis=1) - law_th.central_value

        # sigmas has shape (n_obs, )
        sigmas = np.sqrt(np.var(rep, axis=1))

        delta = bias_diffs / sigmas
        deltas.append(delta.tolist())

    return np.asarray(deltas)


@check_at_least_10_fits
def n_fit_samples(fits):
    """Return a range object where each item is a number of fits to use for
    resampling a multiclosure quantity.

    It is determined by varying n_fits between 10 and number of fits provided by
    user in steps of 5. User must provide at least 10 fits.

    """
    return list(range(10, len(fits) + SAMPLING_INTERVAL, SAMPLING_INTERVAL))


# NOTE: check_multifit_replicas can fill in _internal_max_reps and
# _internal_min_reps if they are None which means by default the values are
# filled but this value can be overridden in specific studies. Both keys
# must be present in signature for the check to work
@check_multifit_replicas
def n_replica_samples(fits_pdf, _internal_max_reps=None, _internal_min_reps=20):
    """Return a range object where each item is a number of replicas to use for
    resampling a multiclosure quantity.

    It is determined by varying n_reps between 20 and number of replicas that each
    provided closure fit has. All provided fits must have the same number of
    replicas and that number must be at least 20.

    The number of replicas used from each fit can be overridden by supplying
    _internal_max_reps.
    """
    return list(
        range(_internal_min_reps, _internal_max_reps + SAMPLING_INTERVAL, SAMPLING_INTERVAL)
    )


class BootstrappedTheoryResult:
    """Proxy class which mimics results.ThPredictionsResult so that
    pre-existing bias/variance actions can be used with bootstrapped replicas
    """

    def __init__(self, data):
        self.error_members = data
        self.central_value = data.mean(axis=1)
        self.rawdata = np.concatenate([self.central_value.reshape(-1, 1), data], axis=-1)


def _bootstrap_multiclosure_fits(
    multiclosure_dataset_loader, rng, n_fit_max, n_fit, n_rep_max, n_rep, use_repeats
):
    """Perform a single bootstrap resample of the multiclosure fits and return
    a proxy of the base internal object used by relevant estimator actions
    with the fits and replicas resampled.

    If use_repeats is False then each fit and replica can only be chosen once
    and there are no repeated samples of either fit or replicas within each fit.

    The various n_fit* and n_rep* choices are for finite size effect studies.
    If you want to perform a simple bootstrap then simply set n_fit and n_fit_max
    to the number of closure fits (len(fits)) and n_rep and n_rep_max to the
    number of replicas in each of the closure tests.

    Returns
    -------

        resampled_multiclosure:
            like multiclosure_dataset_loader but with the fits
            and replicas resampled randomly using np.random.choice.

    See also:

    np.random.choice

    """
    closure_th, *input_tuple = multiclosure_dataset_loader
    fit_boot_index = rng.choice(n_fit_max, size=n_fit, replace=use_repeats)
    fit_boot_th = [closure_th[i] for i in fit_boot_index]
    boot_ths = []
    # construct proxy fits theory predictions
    for fit_th in fit_boot_th:
        rep_boot_index = rng.choice(n_rep_max, size=n_rep, replace=use_repeats)
        boot_ths.append(BootstrappedTheoryResult(fit_th.error_members[:, rep_boot_index]))
    return (boot_ths, *input_tuple)


def bootstrapped_internal_multiclosure_dataset_loader(
    multiclosure_dataset_loader,
    n_fit_max,
    n_fit,
    n_rep_max,
    n_rep,
    n_boot_multiclosure,
    rng_seed_mct_boot,
    use_repeats=True,
):
    """
    Returns a tuple of multiclosure_dataset_loader objects
    each of which is a bootstrap resample of the original dataset

    Parameters
    ----------
    multiclosure_dataset_loader: tuple
        closure fits theory predictions,
        underlying law theory predictions,
        covariance matrix,
        sqrt covariance matrix

    n_fit_max: int
        maximum number of fits, should be smaller or equal to number of
        multiclosure fits

    n_fit: int
        number of fits to draw for each resample

    n_rep_max: int
        maximum number of replicas, should be smaller or equal to number of
        replicas in each fit

    n_rep: int
        number of replicas to draw for each resample

    n_boot_multiclosure: int
        number of bootstrap resamples to perform

    rng_seed_mct_boot: int
        seed for random number generator

    use_repeats: bool, default is True
        whether to allow repeated fits and replicas in each resample

    Returns
    -------
    resampled_multiclosure: tuple of shape (n_boot_multiclosure,)
        tuple of multiclosure_dataset_loader objects each of which
        is a bootstrap resample of the original dataset

    """
    rng = np.random.RandomState(seed=rng_seed_mct_boot)
    return tuple(
        [
            _bootstrap_multiclosure_fits(
                multiclosure_dataset_loader,
                rng=rng,
                n_fit_max=n_fit_max,
                n_fit=n_fit,
                n_rep_max=n_rep_max,
                n_rep=n_rep,
                use_repeats=use_repeats,
            )
            for _ in range(n_boot_multiclosure)
        ]
    )


def bootstrapped_internal_multiclosure_data_loader(
    multiclosure_data_loader,
    n_fit_max,
    n_fit,
    n_rep_max,
    n_rep,
    n_boot_multiclosure,
    rng_seed_mct_boot,
    use_repeats=True,
):
    """Like bootstrapped_internal_multiclosure_dataset_loader except for all data"""
    return bootstrapped_internal_multiclosure_dataset_loader(
        multiclosure_data_loader,
        n_fit_max,
        n_fit,
        n_rep_max,
        n_rep,
        n_boot_multiclosure,
        rng_seed_mct_boot,
        use_repeats,
    )


@check_multifit_replicas
def fits_bootstrap_data_bias_variance(
    multiclosure_data_loader,
    fits,
    _internal_max_reps=None,
    _internal_min_reps=20,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
):
    """Perform bootstrap resample of `fits_data_bias_variance`, returns
    tuple of bias_samples, variance_samples where each element is a 1-D np.array
    of length bootstrap_samples. The elements of the arrays are bootstrap samples
    of bias and variance respectively.

    """
    # seed same rng so we can aggregate results
    rng = np.random.RandomState(seed=boot_seed)
    bias_boot = []
    variance_boot = []
    for _ in range(bootstrap_samples):
        # use all fits. Use all replicas by default. Allow repeats in resample.
        boot_internal_loader = _bootstrap_multiclosure_fits(
            multiclosure_data_loader,
            rng,
            len(fits),
            len(fits),
            _internal_max_reps,
            _internal_max_reps,
            True,
        )
        # explicitly pass n_rep to fits_dataset_bias_variance so it uses
        # full subsample
        bias, variance, _ = expected_dataset_bias_variance(
            fits_dataset_bias_variance(boot_internal_loader, _internal_max_reps, _internal_min_reps)
        )
        bias_boot.append(bias)
        variance_boot.append(variance)
    return np.array(bias_boot), np.array(variance_boot)


def xq2_dataset_map(
    xq2map_with_cuts, multiclosure_dataset_loader, _internal_max_reps=None, _internal_min_reps=20
):
    """
    Load in a dictionary all the specs of a dataset meaning:
    - ds name
    - ds coords
    - standard deviation (in multiclosure)
    - mean (in multiclosure again)
    - (x,Q^2) coords
    """

    commondata = xq2map_with_cuts.commondata
    coords = xq2map_with_cuts[2]
    central_deltas = fits_normed_dataset_central_delta(multiclosure_dataset_loader)
    # std_devs = np.std(central_deltas, axis = 0)
    std_devs = np.sqrt(np.mean(central_deltas**2, axis=0))
    means = np.mean(central_deltas, axis=0)
    xi = dataset_xi(dataset_replica_and_central_diff(multiclosure_dataset_loader, False))

    # for the case of double-hadronic observables we have 2 (x,Q) for each experimental point
    if coords[0].shape[0] != std_devs.shape[0]:
        std_devs = np.concatenate((std_devs, std_devs))
        means = np.concatenate((means, means))
        xi = np.concatenate((xi, xi))
    return {
        'x_coords': coords[0],
        'Q_coords': coords[1],
        'std_devs': std_devs,
        'name': commondata.name,
        'process': commondata.process_type,
        'means': means,
        'xi': xi,
    }


xq2_data_map = collect("xq2_dataset_map", ("data",))


def standard_indicator_function(standard_variable, nsigma=1):
    """
    Calculate the indicator function for a standardised variable.

    Parameters
    ----------
    standard_variable: np.array
        array of variables that have been standardised: (x - mu)/sigma

    nsigma: float
        number of standard deviations to consider

    Returns
    -------
    np.array
        array of ones and zeros. If 1 then the variable is within nsigma standard deviations
        from the mean, otherwise it is 0.
    """
    return np.array(abs(standard_variable) < nsigma, dtype=int)
