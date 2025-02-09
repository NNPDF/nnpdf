"""
Module for bootstrapping multiclosure fits.

"""

import numpy as np
import pandas as pd

from reportengine import collect

import copy
from validphys.closuretest.multiclosure import (
    MulticlosureLoader,
    mean_covmat_multiclosure,
    bias_dataset,
    bias_data,
    regularized_multiclosure_dataset_loader,
    regularized_multiclosure_data_loader,
    normalized_delta_bias_data,
)


class BootstrappedTheoryResult:
    """
    Proxy class which mimics results.ThPredictionsResult so that
    pre-existing bias/variance actions can be used with bootstrapped replicas
    """

    def __init__(self, data):
        self.error_members = data
        self.central_value = data.mean(axis=1)
        self.rawdata = np.concatenate([self.central_value.reshape(-1, 1), data], axis=-1)


def _bootstrap_multiclosure_fits(
    multiclosure_loader: MulticlosureLoader,
    rng: np.random.RandomState,
    n_fit_max: int,
    n_fit: int,
    n_rep_max: int,
    n_rep: int,
    use_repeats: bool,
):
    """
    Perform a single bootstrap resample of the multiclosure fits and return
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
    resampled_multiclosure_loader:
        MulticlosureLoader with bootstrapped closure theories.
    """
    closure_th = multiclosure_loader.closure_theories
    fit_boot_index = rng.choice(n_fit_max, size=n_fit, replace=use_repeats)
    fit_boot_th = [closure_th[i] for i in fit_boot_index]
    boot_ths = []
    # construct proxy fits theory predictions
    for fit_th in fit_boot_th:
        rep_boot_index = rng.choice(n_rep_max, size=n_rep, replace=use_repeats)
        boot_ths.append(BootstrappedTheoryResult(fit_th.error_members[:, rep_boot_index]))

    resampled_multiclosure_loader = copy.deepcopy(multiclosure_loader)
    # replace closure_theories with bootstrapped theories
    resampled_multiclosure_loader.closure_theories = boot_ths

    # Replace the original replica covmat with the new bootstrapped covmat
    covmat_reps_mean = mean_covmat_multiclosure(boot_ths)
    resampled_multiclosure_loader.covmat_reps_mean = covmat_reps_mean

    return resampled_multiclosure_loader


def bootstrapped_multiclosure_dataset_loader(
    multiclosure_dataset_loader: MulticlosureLoader,
    n_fit_max: int,
    n_fit: int,
    n_rep_max: int,
    n_rep: int,
    n_boot_multiclosure: int,
    use_repeats: bool = True,
):
    """
    Returns a tuple of MulticlosureLoader objects
    each of which is a bootstrap resample of the original dataset.

    Parameters
    ----------
    multiclosure_dataset_loader: MulticlosureLoader

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
        tuple of MulticlosureLoader objects each of which
        is a bootstrap resample of the original dataset

    """
    return tuple(
        [
            _bootstrap_multiclosure_fits(
                multiclosure_dataset_loader,
                rng=np.random.RandomState(seed=i),
                n_fit_max=n_fit_max,
                n_fit=n_fit,
                n_rep_max=n_rep_max,
                n_rep=n_rep,
                use_repeats=use_repeats,
            )
            for i in range(n_boot_multiclosure)
        ]
    )


def bootstrapped_multiclosure_data_loader(
    multiclosure_data_loader: MulticlosureLoader,
    n_fit_max: int,
    n_fit: int,
    n_rep_max: int,
    n_rep: int,
    n_boot_multiclosure: int,
    use_repeats: bool = True,
):
    """Like `bootstrapped_multiclosure_dataset_loader` except for all data."""
    return bootstrapped_multiclosure_dataset_loader(
        multiclosure_data_loader,
        n_fit_max,
        n_fit,
        n_rep_max,
        n_rep,
        n_boot_multiclosure,
        use_repeats,
    )


def bootstrapped_regularized_multiclosure_dataset_loader(
    multiclosure_dataset_loader: MulticlosureLoader,
    n_fit_max: int,
    n_fit: int,
    n_rep_max: int,
    n_rep: int,
    n_boot_multiclosure: int,
    use_repeats: bool = True,
    explained_variance_ratio: float = 0.95,
    _internal_max_reps=None,
    _internal_min_reps=20,
) -> tuple:
    """
    Similar to multiclosure.bootstrapped_multiclosure_dataset_loader but returns
    PCA regularised covariance matrix, where the covariance matrix has been computed
    from the replicas of the theory predictions.

    Returns a tuple of RegularizedMulticlosureLoader objects.
    """

    # get bootstrapped multiclosure dataset loader
    bootstrap_mdl = bootstrapped_multiclosure_dataset_loader(
        multiclosure_dataset_loader,
        n_fit_max=n_fit_max,
        n_fit=n_fit,
        n_rep_max=n_rep_max,
        n_rep=n_rep,
        n_boot_multiclosure=n_boot_multiclosure,
        use_repeats=use_repeats,
    )

    # PCA regularise all the bootstrapped internal multiclosure dataset loaders
    bootstrap_mdl_pca = [
        regularized_multiclosure_dataset_loader(
            mdl, explained_variance_ratio, _internal_max_reps, _internal_min_reps
        )
        for mdl in bootstrap_mdl
    ]
    return tuple(bootstrap_mdl_pca)


def bootstrapped_regularized_multiclosure_data_loader(
    multiclosure_data_loader: MulticlosureLoader,
    n_fit_max: int,
    n_fit: int,
    n_rep_max: int,
    n_rep: int,
    n_boot_multiclosure: int,
    use_repeats: bool = True,
    explained_variance_ratio: float = 0.95,
    _internal_max_reps=None,
    _internal_min_reps=20,
) -> tuple:
    """
    Same as `bootstrapped_regularized_multiclosure_dataset_loader` but for all the data.
    """
    # get bootstrapped internal multiclosure dataset loader
    bootstrap_mdl = bootstrapped_multiclosure_data_loader(
        multiclosure_data_loader,
        n_fit_max=n_fit_max,
        n_fit=n_fit,
        n_rep_max=n_rep_max,
        n_rep=n_rep,
        n_boot_multiclosure=n_boot_multiclosure,
        use_repeats=use_repeats,
    )

    # PCA regularise all the bootstrapped internal multiclosure dataset loaders
    bootstrap_mdl_pca = [
        regularized_multiclosure_data_loader(
            mdl, explained_variance_ratio, _internal_max_reps, _internal_min_reps
        )
        for mdl in bootstrap_mdl
    ]
    return tuple(bootstrap_mdl_pca)


def bootstrapped_bias_dataset(bootstrapped_regularized_multiclosure_dataset_loader, dataset):
    """
    Computes Bias for each bootstrap sample.
    Returns a DataFrame with the results.
    """
    boot_bias_var_samples = []
    for i, boot_mdl_pca in enumerate(bootstrapped_regularized_multiclosure_dataset_loader):
        bias, n_comp = bias_dataset(boot_mdl_pca)

        boot_bias_var_samples.append(
            {"bias": np.mean(bias), "n_comp": n_comp, "dataset": str(dataset), "bootstrap_index": i}
        )

    df = pd.DataFrame.from_records(
        boot_bias_var_samples,
        index="bootstrap_index",
        columns=("bootstrap_index", "dataset", "n_comp", "bias"),
    )

    df.columns = ["dataset", "n_comp", "bias"]
    return df


"""
Collect `bootstrapped_bias_dataset` over all datasets.
"""
bootstrapped_bias_datasets = collect("bootstrapped_bias_dataset", ("data",))


def bootstrapped_bias_data(bootstrapped_regularized_multiclosure_data_loader):
    """
    Computes Bias and Variance for each bootstrap sample.
    Returns a DataFrame with the results.
    """
    boot_bias_var_samples = []
    for i, boot_mdl_pca in enumerate(bootstrapped_regularized_multiclosure_data_loader):
        bias, n_comp = bias_data(boot_mdl_pca)
        boot_bias_var_samples.append(
            {
                "bias": np.mean(bias),
                "err bias": np.std(bias),
                "n_comp": n_comp,
                "data": "Full dataset",
                "bootstrap_index": i,
            }
        )

    df = pd.DataFrame.from_records(
        boot_bias_var_samples,
        index="bootstrap_index",
        columns=("bootstrap_index", "dataset", "n_comp", "bias", "err bias"),
    )

    df.columns = ["dataset", "n_comp", "bias", "err bias"]
    return df


def bootstrapped_normalized_delta_bias_data(bootstrapped_regularized_multiclosure_data_loader):
    """
    Compute the normalized deltas for each bootstrap sample.
    Note: delta is the bias in the diagonal basis.

    Parameters
    ----------
    bootstrapped_regularized_multiclosure_data_loader: list
        list of `RegularizedMulticlosureLoader` objects.

    Returns
    -------
    list
    """
    normalised_deltas = []
    for boot_mdl_pca in bootstrapped_regularized_multiclosure_data_loader:
        normalised_deltas.append(normalized_delta_bias_data(boot_mdl_pca))
    return normalised_deltas


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


def bootstrapped_indicator_function_data(bootstrapped_normalized_delta_bias_data, nsigma=1):
    """
    Compute the indicator function for each bootstrap sample.

    Parameters
    ----------
    bootstrapped_normalized_delta_bias_data: list
        list containing the normalized deltas and the number of principal components.
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
    for boot, ndof in bootstrapped_normalized_delta_bias_data:
        indicator_list.append(standard_indicator_function(boot, nsigma))
        ndof_list.append(ndof)
    return indicator_list, np.mean(np.asarray(ndof_list))
