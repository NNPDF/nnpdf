"""
closuretest/multiclosure.py

Module containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit. The actions
in this module are used to produce results which are plotted in
``multiclosure_output.py``

"""

import numpy as np
import dataclasses

from reportengine import collect
import validphys
from validphys.calcutils import calc_chi2
from validphys.checks import check_use_t0
from validphys.closuretest.closure_checks import (
    check_fits_areclosures,
    check_fits_different_filterseed,
    check_fits_underlying_law_match,
    check_multifit_replicas,
    check_t0pdfset_matches_multiclosure_law,
)
from validphys.results import ThPredictionsResult


@dataclasses.dataclass
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


def mean_covmat_multiclosure(closure_theories: list) -> np.array:
    """
    Computes the 'PDF' covariance matrices obtained from each multiclosure
    fit and averages over them.

    Parameters
    ----------
    closure_theories: list
        list of ThPredictionsResult

    Returns
    -------
    covmat_reps_mean
        np.array
    """
    reps = np.asarray([th.error_members for th in closure_theories])
    # compute the covariance matrix of the theory predictions for each fit
    _covmats = np.array([np.cov(rep, rowvar=True, bias=True) for rep in reps])

    # compute the mean covariance matrix
    covmat_reps_mean = np.mean(_covmats, axis=0)

    return covmat_reps_mean


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

    # compute the mean covariance matrix
    covmat_reps_mean = mean_covmat_multiclosure(closure_theories)

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


@dataclasses.dataclass
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
    std_covmat_reps: np.array
        Square root of diagonal entries of the original covariance matrix.
    """

    pc_basis: np.array
    n_comp: int
    reg_covmat_reps_mean: np.array
    sqrt_reg_covmat_reps_mean: np.array
    std_covmat_reps: np.array


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
            std_covmat_reps=np.sqrt(np.diag(covmat_reps_mean)),
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
            std_covmat_reps=np.sqrt(np.diag(covmat_reps_mean)),
        )

    # compute sqrt of pdf covariance matrix (NOTE: the matrix should be diagonal)
    sqrt_reg_covmat_reps_mean = np.diag(np.sqrt(np.diag(reg_covmat_reps_mean)))

    return RegularizedMulticlosureLoader(
        closure_theories=closures_th,
        law_theory=law_th,
        covmat_reps_mean=covmat_reps_mean,
        pc_basis=pc_basis,
        n_comp=n_comp,
        reg_covmat_reps_mean=reg_covmat_reps_mean,
        sqrt_reg_covmat_reps_mean=sqrt_reg_covmat_reps_mean,
        std_covmat_reps=np.sqrt(np.diag(covmat_reps_mean)),
    )


@check_multifit_replicas
def regularized_multiclosure_data_loader(
    multiclosure_data_loader: MulticlosureLoader,
    explained_variance_ratio=0.99,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """
    Similar to multiclosure.regularized_multiclosure_dataset_loader except for all data.
    In this case we regularize the correlation matrix rather than the covariance matrix,
    the reason for this is that different experiments can have different units.

    Parameters
    ----------
    multiclosure_data_loader: MulticlosureLoader
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
            std_covmat_reps=D,
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
    reg_corrmat_reps_mean = pc_basis.T @ _corrmat_mean @ pc_basis

    if n_comp == 1:
        return RegularizedMulticlosureLoader(
            closure_theories=closures_th,
            law_theory=law_th,
            covmat_reps_mean=covmat_reps_mean,
            pc_basis=pc_basis,
            n_comp=1,
            reg_covmat_reps_mean=reg_corrmat_reps_mean,
            sqrt_reg_covmat_reps_mean=np.sqrt(reg_corrmat_reps_mean),
            std_covmat_reps=D,
        )

    # compute sqrt of pdf correlation matrix (NOTE: the matrix should be diagonal)
    sqrt_reg_corrmat_reps_mean = np.diag(np.sqrt(np.diag(reg_corrmat_reps_mean)))

    return RegularizedMulticlosureLoader(
        closure_theories=closures_th,
        law_theory=law_th,
        covmat_reps_mean=covmat_reps_mean,
        pc_basis=pc_basis,
        n_comp=n_comp,
        reg_covmat_reps_mean=reg_corrmat_reps_mean,
        sqrt_reg_covmat_reps_mean=sqrt_reg_corrmat_reps_mean,
        std_covmat_reps=D,
    )


def compute_normalized_bias(
    regularized_multiclosure_loader: RegularizedMulticlosureLoader, corrmat: bool = False
) -> np.array:
    """
    Compute the normalized bias for a RegularizedMulticlosureLoader object.
    If corrmat is True, the bias is computed assuming that RegularizedMulticlosureLoader
    contains the correlation matrix, this is needed when computing the bias for the entire data.

    Parameters
    ----------
    regularized_multiclosure_loader: RegularizedMulticlosureLoader
    corrmat: bool, default is False

    Returns
    -------
    np.array
        Array of shape len(fits) containing the normalized bias for each fit.
    """
    # TODO
    closure_theories = regularized_multiclosure_loader.closure_theories
    law_theory = regularized_multiclosure_loader.law_theory
    n_comp = regularized_multiclosure_loader.n_comp
    pc_basis = regularized_multiclosure_loader.pc_basis
    sqrt_reg_covmat_reps_mean = regularized_multiclosure_loader.sqrt_reg_covmat_reps_mean
    std_covmat_reps = regularized_multiclosure_loader.std_covmat_reps

    reps = np.asarray([th.error_members for th in closure_theories])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - law_theory.central_value[:, np.newaxis]

    if n_comp == 1:
        # TODO: this bit still needs to be tested
        import IPython

        IPython.embed()
        delta_bias = pc_basis * delta_bias
        if corrmat:
            delta_bias /= std_covmat_reps
        biases = (delta_bias / sqrt_reg_covmat_reps_mean) ** 2

    else:
        if corrmat:
            delta_bias = pc_basis.T @ (delta_bias.T / std_covmat_reps).T
        else:
            delta_bias = pc_basis.T @ delta_bias
        biases = calc_chi2(sqrt_reg_covmat_reps_mean, delta_bias)

    return biases


def bias_dataset(regularized_multiclosure_dataset_loader):
    """
    Computes the normalized bias for a RegularizedMulticlosureLoader object
    for a single dataset.

    Parameters
    ----------
    regularized_multiclosure_dataset_loader : RegularizedMulticlosureLoader

    Returns
    -------
    tuple
        bias_fits
        n_comp
    """
    bias_fits = compute_normalized_bias(regularized_multiclosure_dataset_loader, corrmat=False)
    n_comp = regularized_multiclosure_dataset_loader.n_comp
    return bias_fits / n_comp, n_comp


"""
Collects the bias data for all datasets.
"""
bias_datasets = collect("bias_dataset", ("data",))


def bias_data(regularized_multiclosure_data_loader):
    """
    Similar to `bias_dataset` but for all data.
    """
    bias_fits = compute_normalized_bias(regularized_multiclosure_data_loader, corrmat=True)
    n_comp = regularized_multiclosure_data_loader.n_comp
    return bias_fits / n_comp, n_comp


def normalized_delta_bias_data(
    regularized_multiclosure_data_loader: RegularizedMulticlosureLoader,
) -> tuple:
    """
    Compute for all data only the normalized delta after PCA regularization.

    Parameters
    ----------
    regularized_multiclosure_data_loader : tuple
        Tuple containing the results of multiclosure fits after pca regularization

    Returns
    -------
    tuple
        deltas
        n_comp
    """
    # NOTE: function computes delta assuming RegularizedMulticlosureLoader
    # contains the regularized / diagonal correlation matrix.

    pca_loader = regularized_multiclosure_data_loader
    closures_th = pca_loader.closure_theories
    law_th = pca_loader.law_theory
    reg_covmat_reps_mean = pca_loader.reg_covmat_reps_mean
    sqrt_reg_covmat_reps_mean = pca_loader.sqrt_reg_covmat_reps_mean
    pc_basis = pca_loader.pc_basis
    std_covmat_reps = pca_loader.std_covmat_reps
    n_comp = pca_loader.n_comp

    reps = np.asarray([th.error_members for th in closures_th])

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - law_th.central_value[:, np.newaxis]

    # TODO: need to understand the n_comp case
    if n_comp == 1:
        # For full data we regularize the correlation matrix
        delta_bias = pc_basis.T @ (delta_bias / std_covmat_reps)
        std_deviations = sqrt_reg_covmat_reps_mean

    else:
        # delta_bias = eigenvects.T @ (pc_basis.T @ delta_bias)
        # For full data we regularize the correlation matrix
        delta_bias = pc_basis.T @ (delta_bias.T / std_covmat_reps).T
        # reg_covmat_reps_mean should be diagonal
        std_deviations = np.sqrt(np.diag(reg_covmat_reps_mean))[:, None]

    return (delta_bias / std_deviations).flatten(), n_comp


"""
TODO: the code below should be revised by GDC
"""


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


def xq2_dataset_map(
    xq2map_with_cuts, multiclosure_dataset_loader, _internal_max_reps=None, _internal_min_reps=20
):
    """
    TODO: should be described better.

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


"""
TODO
"""
xq2_data_map = collect("xq2_dataset_map", ("data",))
