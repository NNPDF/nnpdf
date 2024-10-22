"""
closuretest/multiclosure.py

Module containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit. The actions
in this module are used to produce results which are plotted in
``multiclosure_output.py``

"""
import numpy as np
import scipy.linalg as la
import scipy.special as special

from reportengine import collect
from validphys.calcutils import calc_chi2
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

from validphys.kinematics import xq2map_with_cuts

# bootstrap seed default
DEFAULT_SEED = 9689372
# stepsize in fits/replicas to use for finite size bootstraps
SAMPLING_INTERVAL = 5


# TODO: deprecate this at some point
@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
@check_use_t0
@check_t0pdfset_matches_multiclosure_law
def internal_multiclosure_dataset_loader(
    dataset, fits_pdf, multiclosure_underlyinglaw, fits, t0_covmat_from_systematics
):
    """Internal function for loading multiple theory predictions for a given
    dataset and a single covariance matrix using underlying law as t0 PDF,
    which is for use with multiclosure statistical estimators. Avoiding memory
    issues from caching the load function of a group of datasets.

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
    fits: list
        list of closure test fits, used to collect ``fits_pdf``

    Returns
    -------

    multiclosure_results: tuple
        a tuple of length 4 containing all necessary dependencies of multiclosure
        statistical estimators in order:

            closure fits theory predictions,
            underlying law theory predictions,
            covariance matrix,
            sqrt covariance matrix

    Notes
    -----
    This function replicates behaviour found elsewhere in validphys, the reason
    for this is that due to the default caching behaviour one can run into
    memory issues when loading the theory predictions for the amount of fits
    typically used in these studies.

    """
    fits_dataset_predictions = [
        ThPredictionsResult.from_convolution(pdf, dataset) for pdf in fits_pdf
    ]
    fits_underlying_predictions = ThPredictionsResult.from_convolution(
        multiclosure_underlyinglaw, dataset
    )

    sqrt_covmat = la.cholesky(t0_covmat_from_systematics, lower=True)
    # TODO: support covmat reg and theory covariance matrix
    # possibly make this a named tuple
    return (
        fits_dataset_predictions,
        fits_underlying_predictions,
        t0_covmat_from_systematics,
        sqrt_covmat,
    )


@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
@check_t0pdfset_matches_multiclosure_law
@check_use_t0
def internal_multiclosure_data_loader(
    data, fits_pdf, multiclosure_underlyinglaw, fits, dataset_inputs_t0_covmat_from_systematics
):
    """Like `internal_multiclosure_dataset_loader` except for all data"""
    return internal_multiclosure_dataset_loader(
        data, fits_pdf, multiclosure_underlyinglaw, fits, dataset_inputs_t0_covmat_from_systematics
    )


@check_multifit_replicas
def fits_dataset_bias_variance(
    internal_multiclosure_dataset_loader,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """For a single dataset, calculate the bias and variance for each fit
    and return tuple (bias, variance, n_data), where bias and variance are
    1-D arrays of length ``len(fits)``.

    For more information on bias see closuretest.bias_dataset and for more information
    on variance see :py:func:`validphys.closuretest.closure_results.variance_dataset`.

    The fits should each have the same underlying law and t0 PDF, but have
    different filterseeds, so that the level 1 shift is different.

    Can control the number of replicas taken from each fit with
    ``_internal_max_reps``.

    """
    closures_th, law_th, _, sqrtcov = internal_multiclosure_dataset_loader
    # The dimentions here are (fit, data point, replica)
    reps = np.asarray([th.error_members[:, :_internal_max_reps] for th in closures_th])
    # take mean across replicas - since we might have changed no. of reps
    centrals = reps.mean(axis=2)
    # place bins on first axis
    diffs = law_th.central_value[:, np.newaxis] - centrals.T
    biases = calc_chi2(sqrtcov, diffs)
    variances = []
    # this seems slow but breaks for datasets with single data point otherwise
    for i in range(reps.shape[0]):
        diffs = reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
        variances.append(np.mean(calc_chi2(sqrtcov, diffs)))
    return biases, np.asarray(variances), len(law_th)

@check_multifit_replicas
def fits_normed_dataset_central_delta(
    internal_multiclosure_dataset_loader,
    _internal_max_reps=None,
    _internal_min_reps=20):
    """
    For each fit calculate the difference between central expectation value and true val. Normalize this
    value by the variance of the differences between replicas and central expectation value (different
    for each fit but expected to vary only a little). Each observable central exp value is 
    expected to be gaussianly distributed around the true value set by the fakepdf
    """
    closures_th, law_th, exp_cov, sqrtcov = internal_multiclosure_dataset_loader
    # The dimentions here are (fit, data point, replica)
    reps = np.asarray([th.error_members[:, :_internal_max_reps] for th in closures_th])
    # One could mask here some reps in order to avoid redundancy of information
    #TODO

    n_data = len(law_th)
    n_fits = np.shape(reps)[0]
    deltas = []
    # There are n_fits pdf_covariances
    # flag to see whether to eliminate dataset
    for i in range(n_fits):
        bias_diffs = np.mean(reps[i], axis = 1) - law_th.central_value
        # bias diffs in the for loop should have dim = (n_obs)
        sigmas = np.sqrt(np.var(reps[i], axis = 1))
        # var diffs has shape (n_obs,reps)
        bias_diffs = bias_diffs/sigmas
        
        deltas.append(bias_diffs.tolist())
        # biases.shape = (n_fits, n_obs_cut/uncut)
        # variances.shape = (n_fits, n_obs_cut/uncut, reps)
    return np.asarray(deltas)


fits_datasets_bias_variance = collect(
    "fits_dataset_bias_variance", ("data",)
)


def expected_dataset_bias_variance(fits_dataset_bias_variance):
    """For a given dataset calculate the expected bias and variance across fits
    then return tuple (expected bias, expected variance, n_data)

    """
    biases, variances, n_data = fits_dataset_bias_variance
    return np.mean(biases), np.mean(variances), n_data


@check_multifit_replicas
def fits_data_bias_variance(
    internal_multiclosure_data_loader,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """Like `fits_dataset_bias_variance` but for all data"""
    return fits_dataset_bias_variance(
        internal_multiclosure_data_loader, _internal_max_reps, _internal_min_reps
    )


def expected_data_bias_variance(fits_data_bias_variance):
    """Like `expected_dataset_bias_variance` except for all data"""
    return expected_dataset_bias_variance(fits_data_bias_variance)


fits_experiments_bias_variance = collect(
    "fits_data_bias_variance", ("group_dataset_inputs_by_experiment",)
)


def fits_total_bias_variance(fits_experiments_bias_variance):
    """Like `fits_dataset_bias_variance` except for all data, assumes there are
    no inter-experiment correlations. That assumption is broken if a theory
    covariance matrix is used.

    """
    bias_total, variance_total, n_total = np.sum(fits_experiments_bias_variance, axis=0)
    return bias_total, variance_total, n_total


datasets_expected_bias_variance = collect("expected_dataset_bias_variance", ("data",))


experiments_expected_bias_variance = collect(
    "expected_data_bias_variance", ("group_dataset_inputs_by_experiment",)
)


def expected_total_bias_variance(fits_total_bias_variance):
    """Like `expected_dataset_bias_variance` except for all data"""
    return expected_dataset_bias_variance(fits_total_bias_variance)


def dataset_replica_and_central_diff(internal_multiclosure_dataset_loader, diagonal_basis=True):
    """For a given dataset calculate sigma, the RMS difference between
    replica predictions and central predictions, and delta, the difference
    between the central prediction and the underlying prediction.

    If ``diagonal_basis`` is ``True`` he differences are calculated in the
    basis which would diagonalise the dataset's covariance matrix. This is the
    default behaviour.

    """
    closures_th, law_th, covmat, _ = internal_multiclosure_dataset_loader
    replicas = np.asarray([th.error_members for th in closures_th])
    centrals = np.mean(replicas, axis=-1)
    underlying = law_th.central_value

    _, e_vec = la.eigh(covmat)

    central_diff = centrals - underlying[np.newaxis, :]
    var_diff_sqrt = centrals[:, :, np.newaxis] - replicas

    if diagonal_basis:
        # project into basis which diagonalises covariance matrix
        var_diff_sqrt = e_vec.T @ var_diff_sqrt.transpose(2, 1, 0)
        central_diff = e_vec.T @ central_diff.T
    else:
        var_diff_sqrt = var_diff_sqrt.transpose(2, 1, 0)
        central_diff = central_diff.T

    var_diff = var_diff_sqrt**2
    sigma = np.sqrt(var_diff.mean(axis=0))  # sigma is always positive
    return sigma, central_diff


def dataset_xi(dataset_replica_and_central_diff):
    """Take sigma and delta for a dataset, where sigma is the RMS difference
    between replica predictions and central predictions, and delta is the
    difference between the central prediction and the underlying prediction.

    Then the indicator function is evaluated elementwise for sigma and delta

        :math:`I_{[-\sigma_j, \sigma_j]}(\delta_j)`

    which is 1 when :math:`|\delta_j| < \sigma_j` and 0 otherwise. Finally, take the
    mean across fits.

    Returns
    -------

        xi_1sigma_i: np.array
            a 1-D array where each element is the value of xi_1sigma for that
            particular eigenvector. We note that the eigenvectors are ordered by
            ascending eigenvalues

    """

    sigma, central_diff = dataset_replica_and_central_diff
    # sigma is always positive
    in_1_sigma = np.array(abs(central_diff) < sigma, dtype=int)
    # mean across fits
    return in_1_sigma.mean(axis=1)


def data_replica_and_central_diff(internal_multiclosure_data_loader, diagonal_basis=True):
    """Like ``dataset_replica_and_central_diff`` but for all data"""
    return dataset_replica_and_central_diff(internal_multiclosure_data_loader, diagonal_basis)


def data_xi(data_replica_and_central_diff):
    """Like dataset_xi but for all data"""
    return dataset_xi(data_replica_and_central_diff)


experiments_xi_measured = collect("data_xi", ("group_dataset_inputs_by_experiment",))
experiments_replica_central_diff = collect(
    "data_replica_and_central_diff", ("group_dataset_inputs_by_experiment",)
)


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
        range(
            _internal_min_reps,
            _internal_max_reps + SAMPLING_INTERVAL,
            SAMPLING_INTERVAL,
        )
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
    internal_multiclosure_dataset_loader,
    rng,
    n_fit_max,
    n_fit,
    n_rep_max,
    n_rep,
    use_repeats,
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
            like internal_multiclosure_dataset_loader but with the fits
            and replicas resampled randomly using np.random.choice.

    See also:

    np.random.choice

    """
    closure_th, *input_tuple = internal_multiclosure_dataset_loader
    fit_boot_index = rng.choice(n_fit_max, size=n_fit, replace=use_repeats)
    fit_boot_th = [closure_th[i] for i in fit_boot_index]
    boot_ths = []
    # construct proxy fits theory predictions
    for fit_th in fit_boot_th:
        rep_boot_index = rng.choice(n_rep_max, size=n_rep, replace=use_repeats)
        boot_ths.append(BootstrappedTheoryResult(fit_th.error_members[:, rep_boot_index]))
    return (boot_ths, *input_tuple)


def bootstrapped_internal_multiclosure_dataset_loader(
    internal_multiclosure_dataset_loader, n_fit_max, n_fit, n_rep_max, n_rep, n_boot_multiclosure, rng_seed_mct_boot, use_repeats=True, 
):
    """
    Returns a tuple of internal_multiclosure_dataset_loader objects
    each of which is a bootstrap resample of the original dataset

    Parameters
    ----------
    internal_multiclosure_dataset_loader: tuple
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
        tuple of internal_multiclosure_dataset_loader objects each of which
        is a bootstrap resample of the original dataset

    """
    rng = np.random.RandomState(seed=rng_seed_mct_boot)
    return tuple([
        _bootstrap_multiclosure_fits(
            internal_multiclosure_dataset_loader,
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


def bias_variance_resampling_dataset(
    internal_multiclosure_dataset_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """For a single dataset, create bootstrap distributions of bias and variance
    varying the number of fits and replicas drawn for each resample. Return two
    3-D arrays with dimensions

        (number of n_rep samples, number of n_fit samples, n_boot)

    filled with resampled bias and variance respectively. The number of bootstrap_samples
    is 100 by default. The number of n_rep samples is determined by varying
    n_rep between 10 and the number of replicas each fit has in intervals of 5.
    This action requires that each fit has the same number of replicas which also
    must be at least 10. The number of n_fit samples is determined analogously to
    the number of n_rep samples, also requiring at least 10 fits.

    Returns
    -------
    resamples: tuple
        tuple of two 3-D arrays with resampled bias and variance respectively for
        each n_rep samples and each n_fit samples

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple datasets then the set of resamples all used corresponding replicas
    and fits.

    """
    # seed same rng so we can aggregate results across datasets
    rng = np.random.RandomState(seed=boot_seed)
    bias_sample = []
    variance_sample = []
    for n_rep_sample in n_replica_samples:
        # results varying n_fit_sample
        fixed_n_rep_bias = []
        fixed_n_rep_variance = []
        for n_fit_sample in n_fit_samples:
            # for each n_fit and n_replica sample store result of each boot resample
            bias_boot = []
            variance_boot = []
            for _ in range(bootstrap_samples):
                boot_internal_loader = _bootstrap_multiclosure_fits(
                    internal_multiclosure_dataset_loader,
                    rng,
                    n_fit_samples[-1],
                    n_fit_sample,
                    n_replica_samples[-1],
                    n_rep_sample,
                    use_repeats,
                )
                # explicitly pass n_rep to fits_dataset_bias_variance so it uses
                # full subsample
                bias, variance, _ = expected_dataset_bias_variance(
                    fits_dataset_bias_variance(boot_internal_loader, n_rep_sample)
                )
                bias_boot.append(bias)
                variance_boot.append(variance)
            fixed_n_rep_bias.append(bias_boot)
            fixed_n_rep_variance.append(variance_boot)
        bias_sample.append(fixed_n_rep_bias)
        variance_sample.append(fixed_n_rep_variance)
    return np.array(bias_sample), np.array(variance_sample)


def bias_variance_resampling_data(
    internal_multiclosure_data_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """Like ratio_n_dependence_dataset except for all data.

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple experiments then the set of resamples all used corresponding
    fits/replicas and can be added together.

    """
    return bias_variance_resampling_dataset(
        internal_multiclosure_data_loader,
        n_fit_samples,
        n_replica_samples,
        bootstrap_samples,
        boot_seed=boot_seed,
        use_repeats=use_repeats,
    )


exps_bias_var_resample = collect(
    "bias_variance_resampling_data", ("group_dataset_inputs_by_experiment",)
)


def bias_variance_resampling_total(exps_bias_var_resample):
    """Sum the bias_variance_resampling_data for all experiments, giving
    the total bias and variance resamples. This relies on the bootstrap seed being
    the same for all experiments, such that the fits/replicas are the same, and
    there being no inter-experiment correlations.

    """
    bias_total, var_total = np.sum(exps_bias_var_resample, axis=0)
    return bias_total, var_total


def xi_resampling_dataset(
    internal_multiclosure_dataset_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """For a single dataset, create bootstrap distributions of xi_1sigma
    varying the number of fits and replicas drawn for each resample. Return a
    4-D array with dimensions

        (number of n_rep samples, number of n_fit samples, n_boot, n_data)

    filled with resampled bias and variance respectively. The number of bootstrap_samples
    is 100 by default. The number of n_rep samples is determined by varying
    n_rep between 10 and the number of replicas each fit has in intervals of 5.
    This action requires that each fit has the same number of replicas which also
    must be at least 10. The number of n_fit samples is determined analogously to
    the number of n_rep samples, also requiring at least 10 fits.

    Returns
    -------
    resamples: array
        4-D array with resampled xi for each n_rep samples and each n_fit samples

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple datasets then the set of resamples all used corresponding replicas.

    """
    # seed same rng so we can aggregate results
    rng = np.random.RandomState(seed=boot_seed)

    xi_1sigma = []
    for n_rep_sample in n_replica_samples:
        # results varying n_fit_sample
        fixed_n_rep_xi_1sigma = []
        for n_fit_sample in n_fit_samples:
            # for each n_fit and n_replica sample store result of each boot resample
            xi_1sigma_boot = []
            for _ in range(bootstrap_samples):
                boot_internal_loader = _bootstrap_multiclosure_fits(
                    internal_multiclosure_dataset_loader,
                    rng,
                    n_fit_samples[-1],
                    n_fit_sample,
                    n_replica_samples[-1],
                    n_rep_sample,
                    use_repeats,
                )
                # append the 1d array for individual eigenvectors
                xi_1sigma_boot.append(
                    dataset_xi(dataset_replica_and_central_diff(boot_internal_loader))
                )
            fixed_n_rep_xi_1sigma.append(xi_1sigma_boot)
        xi_1sigma.append(fixed_n_rep_xi_1sigma)
    return np.array(xi_1sigma)


def xi_resampling_data(
    internal_multiclosure_data_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """Like xi_resampling_dataset except for all data.

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple experiments then the set of resamples all used corresponding replicas
    and can be added together.

    """
    return xi_resampling_dataset(
        internal_multiclosure_data_loader,
        n_fit_samples,
        n_replica_samples,
        bootstrap_samples,
        boot_seed=boot_seed,
        use_repeats=use_repeats,
    )


exps_xi_resample = collect("xi_resampling_data", ("group_dataset_inputs_by_experiment",))


def total_xi_resample(exps_xi_resample):
    """Concatenate the xi for each datapoint for all data"""
    return np.concatenate(exps_xi_resample, axis=-1)


def total_expected_xi_resample(bias_variance_resampling_total):
    """Using the bias and variance resample, return a resample of expected xi
    using the method outlined in
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`.

    The general concept is based on assuming all of the distributions are
    gaussians and using the ratio of bias/variance to predict the corresponding
    integral. To see a more in depth explanation, see
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`.

    """
    bias_total, var_total = bias_variance_resampling_total
    sqrt_bias_var = np.sqrt(bias_total / var_total)
    n_sigma_in_variance = 1 / sqrt_bias_var
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    return special.erf(n_sigma_in_variance / np.sqrt(2))


@check_multifit_replicas
def fits_bootstrap_data_bias_variance(
    internal_multiclosure_data_loader,
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
            internal_multiclosure_data_loader,
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


experiments_bootstrap_bias_variance = collect(
    "fits_bootstrap_data_bias_variance", ("group_dataset_inputs_by_experiment",)
)


def total_bootstrap_ratio(experiments_bootstrap_bias_variance):
    """Calculate the total bootstrap ratio for all data. Leverages the
    fact that the covariance matrix is block diagonal in experiments so

        Total ratio = sum(bias) / sum(variance)

    Which is valid provided there are no inter-experimental correlations.

    Returns
    -------
    bias_var_total: tuple
        tuple of the total bias and variance

    """
    bias_tot, var_tot = np.sum(experiments_bootstrap_bias_variance, axis=0)
    return bias_tot, var_tot


def experiments_bootstrap_ratio(experiments_bootstrap_bias_variance, total_bootstrap_ratio):
    """Returns a bootstrap resampling of the ratio of bias/variance for
    each experiment and total. Total is calculated as sum(bias)/sum(variance)
    where each sum refers to the sum across experiments.

    Returns
    -------

    ratios_resampled: list
        list of bootstrap samples of ratio of bias/var, length of list is
        len(experiments) + 1 because the final element is the total ratio
        resampled.

    """
    ratios = [bias / var for bias, var in experiments_bootstrap_bias_variance]
    bias_tot, var_tot = total_bootstrap_ratio
    ratios.append(bias_tot / var_tot)
    return ratios


def experiments_bootstrap_sqrt_ratio(experiments_bootstrap_ratio):
    """Square root of experiments_bootstrap_ratio"""
    return np.sqrt(experiments_bootstrap_ratio)


def experiments_bootstrap_expected_xi(experiments_bootstrap_sqrt_ratio):
    """Calculate a bootstrap resampling of the expected xi from
    ``experiments_bootstrap_sqrt_ratio``, using the same formula as
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`.

    """
    n_sigma_in_variance = 1 / experiments_bootstrap_sqrt_ratio
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    estimated_integral = special.erf(n_sigma_in_variance / np.sqrt(2))
    return estimated_integral


groups_bootstrap_bias_variance = collect(
    "fits_bootstrap_data_bias_variance", ("group_dataset_inputs_by_metadata",)
)


def groups_bootstrap_ratio(groups_bootstrap_bias_variance, total_bootstrap_ratio):
    """Like :py:func:`experiments_bootstrap_ratio` but for metadata groups."""
    return experiments_bootstrap_ratio(groups_bootstrap_bias_variance, total_bootstrap_ratio)


def groups_bootstrap_sqrt_ratio(groups_bootstrap_ratio):
    """Like :py:func:`experiments_bootstrap_sqrt_ratio` but for metadata groups."""
    return experiments_bootstrap_sqrt_ratio(groups_bootstrap_ratio)


def groups_bootstrap_expected_xi(groups_bootstrap_sqrt_ratio):
    """Like :py:func:`experiments_bootstrap_expected_xi` but for metadata groups."""
    return experiments_bootstrap_expected_xi(groups_bootstrap_sqrt_ratio)


@check_multifit_replicas
def fits_bootstrap_data_xi(
    internal_multiclosure_data_loader,
    fits,
    _internal_max_reps=None,
    _internal_min_reps=20,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
):
    """Perform bootstrap resample of ``data_xi``, returns a list
    where each element is an independent resampling of ``data_xi``.

    For more information on bootstrapping see _bootstrap_multiclosure_fits.
    For more information on xi see dataset_xi.

    """
    # seed same rng so we can aggregate results
    rng = np.random.RandomState(seed=boot_seed)

    xi_1sigma_boot = []
    for _ in range(bootstrap_samples):
        # use all fits. Use all replicas by default. Allow repeats in resample.
        boot_internal_loader = _bootstrap_multiclosure_fits(
            internal_multiclosure_data_loader,
            rng,
            len(fits),
            len(fits),
            _internal_max_reps,
            _internal_max_reps,
            True,
        )
        xi_1sigma_boot.append(dataset_xi(dataset_replica_and_central_diff(boot_internal_loader)))
    return xi_1sigma_boot


experiments_bootstrap_xi = collect(
    "fits_bootstrap_data_xi", ("group_dataset_inputs_by_experiment",)
)


def total_bootstrap_xi(experiments_bootstrap_xi):
    """Given the bootstrap samples of xi_1sigma for all experiments,
    concatenate the result to get xi_1sigma for all data points in a single
    array

    """
    return np.concatenate(experiments_bootstrap_xi, axis=1)


groups_bootstrap_xi = collect("fits_bootstrap_data_xi", ("group_dataset_inputs_by_metadata",))


def dataset_fits_bias_replicas_variance_samples(
    internal_multiclosure_dataset_loader,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """For a single dataset, calculate the samples of chi2-quantities which
    are used to calculate the bias and variance for each fit. The output of this
    function is similar to :py:func:`fits_dataset_bias_variance` except that
    the mean is not taken across replicas when calculating the mean squared
    difference between replica predictions and central predictions and instead
    the results are concatenated. The mean of this array would be the expected
    value of the variance across fits.

    Return tuple (fits_bias, fits_replica_variance, n_data), where fits_bias is
    1-D array of length N_fits and fits_replica_variance is 1-D array length
    N_fits * N_replicas.

    For more information on bias see closuretest.bias_dataset and for more information
    on variance see :py:func:`validphys.closuretest.closure_results.variance_dataset`.

    The fits should each have the same underlying law and t0 PDF, but have
    different filterseeds, so that the level 1 shift is different.

    Can control the number of replicas taken from each fit with
    ``_internal_max_reps``.

    """
    closures_th, law_th, _, sqrtcov = internal_multiclosure_dataset_loader
    # The dimentions here are (fit, data point, replica)
    reps = np.asarray([th.error_members[:, :_internal_max_reps] for th in closures_th])
    # take mean across replicas - since we might have changed no. of reps
    centrals = reps.mean(axis=2)
    # place bins on first axis
    diffs = law_th.central_value[:, np.newaxis] - centrals.T
    biases = calc_chi2(sqrtcov, diffs)
    variances = []
    # this seems slow but breaks for datasets with single data point otherwise
    for i in range(reps.shape[0]):
        diffs = reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
        variances.append(calc_chi2(sqrtcov, diffs))
    return biases, np.concatenate(variances), len(law_th)


def dataset_inputs_fits_bias_replicas_variance_samples(
    internal_multiclosure_data_loader,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    return dataset_fits_bias_replicas_variance_samples(
        internal_multiclosure_data_loader,
        _internal_max_reps=None,
        _internal_min_reps=20,
    )


experiments_fits_bias_replicas_variance_samples = collect(
    "dataset_inputs_fits_bias_replicas_variance_samples", ("group_dataset_inputs_by_experiment",)
)


def xq2_dataset_map(commondata, cuts,internal_multiclosure_dataset_loader,
                        _internal_max_reps=None,
                        _internal_min_reps=20):
    """
    Load in a dictionary all the specs of a dataset meaning:
    - ds name
    - ds coords
    - standard deviation (in multiclosure)
    - mean (in multiclosure again)
    - (x,Q^2) coords
    """
    xq2_map_obj = xq2map_with_cuts(commondata, cuts)
    coords = xq2_map_obj[2]
    central_deltas = fits_normed_dataset_central_delta(internal_multiclosure_dataset_loader)
    std_devs = np.std(central_deltas, axis = 0)
    means = np.mean(central_deltas, axis = 0)
    xi = dataset_xi(dataset_replica_and_central_diff(internal_multiclosure_dataset_loader,False))

    # for case of DY observables we have 2 (x,Q) for each experimental point
    if coords[0].shape[0] != std_devs.shape[0]:
        std_devs = np.concatenate((std_devs,std_devs))
        means = np.concatenate((means,means))
        xi = np.concatenate((xi,xi))
    return {'x_coords':coords[0], 'Q_coords':coords[1],'std_devs':std_devs,'name':commondata.name,'process':commondata.process_type, 'means': means, 'xi': xi}


xq2_data_map = collect("xq2_dataset_map",("data",))