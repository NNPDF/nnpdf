"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing

from validphys import covmats
from validphys.calcutils import calc_chi2
from validphys.results import ThPredictionsResult

from reportengine import collect


""" To load several multiclosure fits. Useful for inconsistent closure test analysis """
multi_dataset_loader = collect("internal_multiclosure_dataset_loader", ("dataspecs",))

multi_dataset_fits_bias_replicas_variance_samples = collect(
    "dataset_fits_bias_replicas_variance_samples", ("dataspecs",)
)

multi_fits_bootstrap_dataset_bias_variance = collect(
    "fits_bootstrap_dataset_bias_variance", ("dataspecs",)
)

multi_bias_variance_resampling_dataset = collect("bias_variance_resampling_dataset", ("dataspecs",))

multi_dataset_fits_bias_variance_samples_pca = collect(
    "dataset_fits_bias_variance_samples_pca", ("dataspecs",)
)


def principal_components_dataset(dataset, fits_pdf, variancepdf, explained_variance_ratio=0.93):
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
        2D tuple:
        - matrix of the principal components (PCs) of shape (N_pc, N_dat)
        - reduced feature matrix, i.e., feature matrix projected onto PCs of shape (N_pc, N_rep)

    """
    # fits_dataset_predictions = [
    #     ThPredictionsResult.from_convolution(pdf, dataset) for pdf in fits_pdf
    # ]

    # dimensions here are (Nfits, Ndat, Nrep)
    # reps = np.asarray([th.error_members for th in fits_dataset_predictions])

    # reshape so as to get PCs from all the samples
    # reps = reps.reshape(reps.shape[1],-1)

    # get replicas from variance fit, used to estimate variance
    reps = ThPredictionsResult.from_convolution(variancepdf, dataset).error_members

    # rescale feature matrix
    reps_scaled = reps  # preprocessing.scale(reps)

    # choose number of principal components (PCs) based on explained variance ratio
    n_comp = 1
    for _ in range(reps.shape[0]):
        pca = PCA(n_comp).fit(reps_scaled.T)
        if np.sum(pca.explained_variance_ratio_) >= explained_variance_ratio:
            break
        n_comp += 1

    # project feature matrix onto PCs
    pc_reps = pca.components_ @ reps

    return pca.components_, pc_reps, n_comp


def principal_components_bias_variance_dataset(
    internal_multiclosure_dataset_loader, principal_components_dataset
):
    """
    TODO
    """

    closures_th, law_th, _, _ = internal_multiclosure_dataset_loader

    reps = np.asarray([th.error_members for th in closures_th])

    pc_basis, pc_reps, n_comp = principal_components_dataset

    # estimate (PC) pdf covariance matrix (from replicas), shape is (Npc, Npc)
    covmat_pdf = np.cov(pc_reps)
    sqrt_covmat_pdf = covmats.sqrt_covmat(covmat_pdf)

    # compute bias diff and project it onto space spanned by PCs
    delta_bias = reps.mean(axis=2).T - law_th.central_value[:, np.newaxis]
    # shape here is (Npc, Nfits)
    delta_bias = pc_basis @ delta_bias

    # compute biases, shape of biases is (Nfits)
    biases = calc_chi2(sqrt_covmat_pdf, delta_bias)

    variances = []
    for i in range(reps.shape[0]):
        diffs = pc_basis @ (reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True))
        variances.append(np.mean(calc_chi2(sqrt_covmat_pdf, diffs)))

    return biases, np.asarray(variances), n_comp

principal_components_bias_variance_datasets = collect(
    "principal_components_bias_variance_dataset", ("data",)
)

def compute_num_components(covariance_matrix, threshold=0.99):
    """
    Compute the number of principal components to keep based on a desired explained variance threshold.

    Parameters
    ----------
    covariance_matrix : np.ndarray, 2-D array

    threshold : (float): Desired explained variance threshold (between 0 and 1).

    Returns
    -------
    int
        num_components: Number of principal components to keep.

    """
    eig_val, eig_vec = np.linalg.eig(covariance_matrix)
    idx = eig_val.argsort()[::-1]
    eig_val = eig_val[idx]
    eig_vec = eig_vec[:, idx]

    cumulative_sum = np.cumsum(eig_val)
    total_sum = np.sum(eig_val)
    num_components = np.argmin(np.abs(cumulative_sum / total_sum - threshold))

    return num_components


def pca_covmat(X, num_components):
    """
    given data X of shape (n,p), reduce its dimension to
    (n,num_components) and return the covariance matrix
    of the reduced data matrix.

    Parameters
    ----------

    Returns
    -------
    """
    pca = PCA(num_components)
    X_reduced = pca.fit_transform(X)
    covariance = np.dot(X_reduced.T, X_reduced) / (X_reduced.shape[0] - 1)
    return covariance


def calc_chi2_pca(pdf_cov, diff, num_components):
    """
    Computes the chi2 between pdf_cov and diff by first projecting
    them into num_components PCs.

    Parameters
    ----------
    pdf_cov: np.ndarray

    diff: np.ndarray

    num_components: int

    Returns
    -------
    float or np.ndarray depending on input

    """
    # Diagonalize the matrix
    L, W = np.linalg.eig(pdf_cov)

    # Sort the eigenvalues and eigenvectors from largest to smallest
    idx = L.argsort()[::-1]
    L = L[idx]
    W = W[:, idx]

    # Keep only the n_components largest eigenvectors
    Wtilde = W[:, :num_components]

    # Transform initial covariance matrix
    pdf_cov_pca = np.einsum("ij,jk->ik", np.einsum("ij,ik->jk", Wtilde, pdf_cov), Wtilde).real

    # transform data
    diff_pca = diff.T @ Wtilde

    try:
        sqrt_pdf_cov_pca = covmats.sqrt_covmat(pdf_cov_pca)
    except:
        raise ValueError(
            f"PCA Covariance Matrix cannot be Cholesky decomposed, perhaps less than {num_components} PC should be kept"
        )

    return np.real(calc_chi2(sqrt_pdf_cov_pca, diff_pca.T))


def dataset_fits_bias_variance_samples_pca(internal_multiclosure_dataset_loader, threshold=0.99):
    """
    similar to `dataset_fits_bias_replicas_variance_samples`.

    Returns
    -------
    tuple
        3D tuple:
        - biases: 1-D array of shape (Nfits,)
        - variances: 1-D array of shape (Nfits, )
        - n_eig: number of eigenvalues kept in PCA, i.e.,
          ndata in the new, lower dimensional, space.
    Note that we return Nfits values of the variance so that computing the
    Bias - Variance ratio is straightforward.
    """
    closures_th, law_th, _, _ = internal_multiclosure_dataset_loader

    # The dimensions here are (fit, data point, replica)
    reps = np.asarray([th.error_members for th in closures_th])

    # take mean across replicas - since we might have changed no. of reps
    centrals = reps.mean(axis=2)

    # compute the PDF covariance matrix of the central samples
    if centrals.shape[0] <= 1:
        raise ValueError(f"Need more than one fit to compute the 'Bias' PDF covariance Matrix")

    pdf_cov_bias = np.cov(centrals.T)

    # find number of (ordered) eigenvalues that explain 99% of the total variance (total sum of eigenvalues)
    n_eig = compute_num_components(pdf_cov_bias, threshold)

    # compute bias from PCs
    diffs_bias = centrals.T - law_th.central_value[:, np.newaxis]
    biases = calc_chi2_pca(pdf_cov_bias, diffs_bias, n_eig)

    # compute variances from PCs
    variances = []

    # loop over fits to compute variances
    for i in range(reps.shape[0]):
        diffs_var = reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
        pdf_cov_var = np.cov(reps[i, :, :])

        variances.append(np.mean(calc_chi2_pca(pdf_cov_var, diffs_var, n_eig)))

    return biases, np.asarray(variances), n_eig


def expected_dataset_fits_bias_variance_samples_pca(dataset_fits_bias_variance_samples_pca):
    """
    returns the bias and variance expected values as well as the number of
    principal components.
    """
    biases, variances, n_eig = dataset_fits_bias_variance_samples_pca
    return np.mean(biases), np.mean(variances), n_eig


expected_datasets_fits_bias_variance_samples_pca = collect(
    "expected_dataset_fits_bias_variance_samples_pca", ("data",)
)


def dataset_fits_ratio_bias_variance_samples_pca(dataset_fits_bias_variance_samples_pca):
    """ """
    biases, variances, _ = dataset_fits_bias_variance_samples_pca
    sqrt_ratios = np.sqrt(biases / variances)
    return sqrt_ratios


def dataset_fits_gaussian_parameters(internal_multiclosure_dataset_loader, threshold=0.99):
    """
    returns parameters of multi gaussian distribution of replicas
    and central replicas
    """
    closures_th, law_th, _, _ = internal_multiclosure_dataset_loader

    # The dimensions here are (fit, data point, replica)
    reps = np.asarray([th.error_members for th in closures_th])

    # take mean across replicas - since we might have changed no. of reps
    centrals = reps.mean(axis=2)

    centrals_covmat = np.cov(centrals.T)
    centrals_covmat = pca_covmat(
        centrals, num_components=compute_num_components(centrals_covmat, threshold)
    )
    mean_centrals = np.mean(centrals, axis=0)

    replicas_covmat = 0
    for i in range(reps.shape[0]):
        replicas_covmat = np.cov(reps[i, :, :])
        replicas_covmat += pca_covmat(
            reps[i, :, :].T, num_components=compute_num_components(replicas_covmat, threshold)
        )
    replicas_covmat /= reps.shape[0]
    mean_replicas = np.mean(reps.reshape(reps.shape[1], -1), axis=1)

    return mean_centrals, centrals_covmat, mean_replicas, replicas_covmat
