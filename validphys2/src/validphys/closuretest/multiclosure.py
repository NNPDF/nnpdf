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

from validphys.results import ThPredictionsResult
from validphys.calcutils import calc_chi2
from validphys.core import DataSetSpec
from validphys.closuretest.closure_checks import (
    check_at_least_10_fits,
    check_multifit_replicas,
    check_fits_underlying_law_match,
    check_fits_areclosures,
    check_fits_different_filterseed,
)

# bootstrap seed default
DEFAULT_SEED = 9689372
# stepsize in fits/replicas to use for finite size bootstraps
SAMPLING_INTERVAL = 5

# TODO: deprecate this at some point
# TODO: update the docstring wording here when ``data`` keyword gets merged.
@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
def internal_multiclosure_dataset_loader(
    dataset, fits_pdf, multiclosure_underlyinglaw, fits
):
    """Internal function for loading multiple theory predictions for a given
    experiment and a single covariance matrix using underlying law as t0 PDF,
    which is for use with multiclosure statistical estimators. Avoiding memory
    issues from caching experiment load.

    Parameters
    ----------
    dataset: DatasetSpec-like
        dataset for which the theory predictions and t0 covariance matrix
        will be loaded. Note that due to the structure of `validphys` this
        function can be overloaded to accept an ExperimentSpec.
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
    if isinstance(dataset, DataSetSpec):
        data = dataset.load()  # just use internal loader
    else:
        # don't cache result
        data = dataset.load.__wrapped__(dataset)

    fits_dataset_predictions = [
        ThPredictionsResult.from_convolution(pdf, dataset, loaded_data=data)
        for pdf in fits_pdf
    ]
    fits_underlying_predictions = ThPredictionsResult.from_convolution(
        multiclosure_underlyinglaw, dataset, loaded_data=data
    )

    # copy data to make t0 cov
    loaded_data = type(data)(data)
    loaded_data.SetT0(multiclosure_underlyinglaw.load_t0())
    covmat = loaded_data.get_covmat()
    sqrt_covmat = la.cholesky(covmat, lower=True)
    # TODO: support covmat reg and theory covariance matrix
    # possibly make this a named tuple
    return (fits_dataset_predictions, fits_underlying_predictions, covmat, sqrt_covmat)


@check_fits_underlying_law_match
@check_fits_areclosures
@check_fits_different_filterseed
def internal_multiclosure_experiment_loader(
    experiment, fits_pdf, multiclosure_underlyinglaw, fits
):
    """Like `internal_multiclosure_dataset_loader` except for an experiment"""
    return internal_multiclosure_dataset_loader(
        experiment, fits_pdf, multiclosure_underlyinglaw, fits
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
    reps = np.asarray([th._rawdata[:, :_internal_max_reps] for th in closures_th])
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


def expected_dataset_bias_variance(fits_dataset_bias_variance):
    """For a given dataset calculate the expected bias and variance across fits
    then return tuple (expected bias, expected variance, n_data)

    """
    biases, variances, n_data = fits_dataset_bias_variance
    return np.mean(biases), np.mean(variances), n_data


@check_multifit_replicas
def fits_experiment_bias_variance(
    internal_multiclosure_experiment_loader,
    _internal_max_reps=None,
    _internal_min_reps=20,
):
    """Like `fits_dataset_bias_variance` but for an experiment"""
    return fits_dataset_bias_variance(
        internal_multiclosure_experiment_loader, _internal_max_reps, _internal_min_reps
    )


def expected_experiment_bias_variance(fits_experiment_bias_variance):
    """Like `expected_dataset_bias_variance` except for an experiment"""
    return expected_dataset_bias_variance(fits_experiment_bias_variance)


fits_experiments_bias_variance = collect(
    "fits_experiment_bias_variance", ("experiments",)
)

# TODO: get rid of this with data keyword merge
def fits_total_bias_variance(fits_experiments_bias_variance):
    """Like `fits_dataset_bias_variance` except for all data"""
    bias_total, variance_total, n_total = np.sum(fits_experiments_bias_variance, axis=0)
    return bias_total, variance_total, n_total


datasets_expected_bias_variance = collect(
    "expected_dataset_bias_variance", ("experiments", "experiment")
)



experiments_expected_bias = collect("expected_bias_experiment", ("experiments",))
experiments_expected_variance = collect(
    "expected_variance_experiment", ("experiments",)
)
experiments_expected_bias_variance = collect(
    "expected_experiment_bias_variance", ("experiments",)
)


def expected_total_bias_variance(fits_total_bias_variance):
    """Like `expected_dataset_bias_variance` except for all data"""
    return expected_dataset_bias_variance(fits_total_bias_variance)



def dataset_xi(internal_multiclosure_dataset_loader):
    """For a given dataset calculate sigma, the RMS difference between
    replica predictions and central predictions, and delta, the difference
    between the central prediction and the underlying prediction.

    The differences are calculated in the basis which would diagonalise the
    dataset's covariance matrix.

    Then the indicator function is evaluated elementwise for sigma and delta

        :math:`I_{[-\sigma_j, \sigma_j]}(\delta_j)`

    which is 1 when :math:`|\delta_j| < \sigma_j` and 0 otherwise. Finally, take the
    mean across fits.

    Returns
    -------

        xi_1sigma_i: np.array
            a 1-D array where each element is the value of xi_1sigma for that
            particular direction. We note that the directions are ordered by
            ascending eigenvalues

    """
    closures_th, law_th, covmat, _ = internal_multiclosure_dataset_loader
    replicas = np.asarray([th._rawdata for th in closures_th])
    centrals = np.mean(replicas, axis=-1)
    underlying = law_th.central_value

    _, e_vec = la.eigh(covmat)

    central_diff = centrals - underlying[np.newaxis, :]
    var_diff_sqrt = centrals[:, :, np.newaxis] - replicas

    # project into basis which diagonalises covariance matrix
    var_diff_sqrt = e_vec.T @ var_diff_sqrt.transpose(2, 1, 0)
    central_diff = e_vec.T @ central_diff.T

    var_diff = var_diff_sqrt ** 2
    sigma = np.sqrt(var_diff.mean(axis=0))  # sigma is always positive
    in_1_sigma = np.array(abs(central_diff) < sigma, dtype=int)
    # mean across fits
    return in_1_sigma.mean(axis=1)


def experiment_xi(internal_multiclosure_experiment_loader):
    """Like dataset_xi but for whole experiment"""
    return dataset_xi(internal_multiclosure_experiment_loader)


experiments_xi_measured = collect("experiment_xi", ("experiments",))




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
        self._rawdata = data
        self.central_value = data.mean(axis=1)


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
        boot_ths.append(BootstrappedTheoryResult(fit_th._rawdata[:, rep_boot_index]))
    return (boot_ths, *input_tuple)


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


def bias_variance_resampling_experiment(
    internal_multiclosure_experiment_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """Like ratio_n_dependence_dataset except for an experiment.

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple experiments then the set of resamples all used corresponding
    fits/replicas and can be added together.

    """
    return bias_variance_resampling_dataset(
        internal_multiclosure_experiment_loader,
        n_fit_samples,
        n_replica_samples,
        bootstrap_samples,
        boot_seed=boot_seed,
        use_repeats=use_repeats,
    )


exps_bias_var_resample = collect(
    "bias_variance_resampling_experiment", ("experiments",)
)


def bias_variance_resampling_total(exps_bias_var_resample):
    """Sum the bias_variance_resampling_experiment for all experiments, giving
    the total bias and variance resamples. This relies on the bootstrap seed being
    the same for all experiments such that the fits/replicas are the same.

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
                # append the 1d array for individual directions
                xi_1sigma_boot.append(dataset_xi(boot_internal_loader))
            fixed_n_rep_xi_1sigma.append(xi_1sigma_boot)
        xi_1sigma.append(fixed_n_rep_xi_1sigma)
    return np.array(xi_1sigma)


def xi_resampling_experiment(
    internal_multiclosure_experiment_loader,
    n_fit_samples,
    n_replica_samples,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
    use_repeats=True,
):
    """Like xi_resampling_dataset except for an experiment.

    Notes
    -----
    The bootstrap samples are seeded in this function. If this action is collected
    over multiple experiments then the set of resamples all used corresponding replicas
    and can be added together.

    """
    return xi_resampling_dataset(
        internal_multiclosure_experiment_loader,
        n_fit_samples,
        n_replica_samples,
        bootstrap_samples,
        boot_seed=boot_seed,
        use_repeats=use_repeats,
    )




exps_xi_resample = collect("xi_resampling_experiment", ("experiments",))


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
def fits_bootstrap_experiment_bias_variance(
    internal_multiclosure_experiment_loader,
    fits,
    _internal_max_reps=None,
    _internal_min_reps=20,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
):
    """Perform bootstrap resample of `fits_experiment_bias_variance`, returns
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
            internal_multiclosure_experiment_loader,
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
            fits_dataset_bias_variance(
                boot_internal_loader, _internal_max_reps, _internal_min_reps
            )
        )
        bias_boot.append(bias)
        variance_boot.append(variance)
    return np.array(bias_boot), np.array(variance_boot)


experiments_bootstrap_bias_variance = collect(
    "fits_bootstrap_experiment_bias_variance", ("experiments",)
)


def experiments_bootstrap_ratio(experiments_bootstrap_bias_variance):
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
    bias_tot, var_tot = experiments_bootstrap_bias_variance[0]
    # add first ratio to list
    ratios = [bias_tot / var_tot]
    for bias, var in experiments_bootstrap_bias_variance[1:]:
        bias_tot += bias
        var_tot += var
        ratios.append(bias / var)
    ratios.append(bias_tot / var_tot)
    return ratios


def experiments_bootstrap_sqrt_ratio(experiments_bootstrap_ratio):
    """Square root of experiments_bootstrap_ratio"""
    return np.sqrt(experiments_bootstrap_ratio)



def experiments_bootstrap_expected_xi(experiments_bootstrap_sqrt_ratio):
    """Calculate a bootstrap resampling of the expected xi from
    `experiments_bootstrap_sqrt_ratio`, using the same formula as
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`.

    """
    n_sigma_in_variance = 1 / experiments_bootstrap_sqrt_ratio
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    estimated_integral = special.erf(n_sigma_in_variance / np.sqrt(2))
    return estimated_integral




@check_multifit_replicas
def fits_bootstrap_experiment_xi(
    internal_multiclosure_experiment_loader,
    fits,
    _internal_max_reps=None,
    _internal_min_reps=20,
    bootstrap_samples=100,
    boot_seed=DEFAULT_SEED,
):
    """Perform bootstrap resample of `experiment_xi`, returns a list
    where each element is an independent resampling of experiment_xi.

    For more information on bootstrapping see _bootstrap_multiclosure_fits.
    For more information on xi see dataset_xi.

    """
    # seed same rng so we can aggregate results
    rng = np.random.RandomState(seed=boot_seed)

    xi_1sigma_boot = []
    for _ in range(bootstrap_samples):
        # use all fits. Use all replicas by default. Allow repeats in resample.
        boot_internal_loader = _bootstrap_multiclosure_fits(
            internal_multiclosure_experiment_loader,
            rng,
            len(fits),
            len(fits),
            _internal_max_reps,
            _internal_max_reps,
            True,
        )
        xi_1sigma_boot.append(dataset_xi(boot_internal_loader))
    return xi_1sigma_boot


experiments_bootstrap_xi = collect("fits_bootstrap_experiment_xi", ("experiments",))

def total_bootstrap_xi(experiments_bootstrap_xi):
    """Given the bootstrap samples of xi_1sigma for all experiments,
    concatenate the result to get xi_1sigma for all data points in a single
    array

    """
    return np.concatenate(experiments_bootstrap_xi, axis=1)
