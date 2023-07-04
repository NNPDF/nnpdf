"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np

from validphys import covmats
from validphys.calcutils import calc_chi2

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


def dataset_fits_bias_variance_samples_pca(internal_multiclosure_dataset_loader, threshold=0.999):
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
    if centrals.shape[0] <=1:
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
