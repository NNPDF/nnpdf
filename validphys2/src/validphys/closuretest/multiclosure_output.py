"""
multiclosure_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for multiclosure estimators in the space of
data.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.special as special

from reportengine.figure import figure
from reportengine.table import table

from validphys.results import each_dataset


@figure
def plot_dataset_fits_bias_variance(fits_dataset_bias_variance, dataset):
    """For a set of closure fits, calculate the bias and variance across fits
    and then plot scatter points so we can see the distribution of each quantity
    with fits. The spread of the variance across fits is assumed to be small
    compared to the spread of the biases, deviation from this assumption could
    suggest that there are finite size effects due to too few replicas.

    """
    biases, variances, _ = fits_dataset_bias_variance
    fig, ax = plt.subplots()
    ax.plot(biases, "*", label=f"bias, std. dev. = {np.std(biases):.2f}")
    ax.axhline(
        np.mean(biases), label=f"bias, mean = {np.mean(biases):.2f}", linestyle="-"
    )
    ax.plot(variances, ".", label=f"variance, std. dev. = {np.std(variances):.2f}")
    ax.axhline(
        np.mean(variances),
        label=f"variance, mean = {np.mean(variances):.2f}",
        linestyle=":",
    )
    ax.set_title(f"Bias and variance for {dataset} for each fit (unnormalised)")
    ax.set_xlabel("fit index")
    ax.legend()
    return fig


@figure
def plot_experiment_fits_bias_variance(fits_experiment_bias_variance, experiment):
    """Like `plot_dataset_fits_bias_variance` but for an experiment.

    """
    return plot_dataset_fits_bias_variance(fits_experiment_bias_variance, experiment)


@figure
def plot_total_fits_bias_variance(fits_total_bias_variance):
    """Like `plot_dataset_fits_bias_variance` but for the total bias/variance
    for all data.

    """
    return plot_dataset_fits_bias_variance(fits_total_bias_variance, "all data")


@table
def datasets_bias_variance_ratio(datasets_expected_bias_variance, each_dataset):
    """For each dataset calculate the expected bias and expected variance
    across fits, then calculate the ratio

        ratio = expected bias / expected variance

    and tabulate the results.

    This gives an idea of how faithful uncertainties are for a set of
    datasets.

    Notes
    -----

    If uncertainties are faithfully estimated then we would expect to see
    ratio = 1. We should note that the ratio is a squared quantity and
    sqrt(ratio) is more appropriate for seeing how much uncertainties are
    over or underestimated. An over-estimate of uncertainty leads to
    sqrt(ratio) < 1, similarly an under-estimate of uncertainty leads to
    sqrt(ratio) > 1.

    """
    records = []
    for ds, (bias, var, ndata) in zip(each_dataset, datasets_expected_bias_variance):
        records.append(dict(dataset=str(ds), ndata=ndata, ratio=bias / var))
    df = pd.DataFrame.from_records(
        records, index="dataset", columns=("dataset", "ndata", "ratio")
    )
    df.columns = ["ndata", "bias/variance"]
    return df


@table
def experiments_bias_variance_ratio(
    experiments_expected_bias_variance, experiments, expected_total_bias_variance
):
    """Like datasets_bias_variance_ratio except for each experiment. Also
    calculate and tabulate

        (total expected bias) / (total expected variance)

    where the total refers to summing over all experiments.

    """
    # don't reinvent wheel
    df_in = datasets_bias_variance_ratio(
        experiments_expected_bias_variance, experiments
    )

    bias_tot, var_tot, ntotal = expected_total_bias_variance

    tot_df = pd.DataFrame(
        [[ntotal, bias_tot / var_tot]], index=["Total"], columns=df_in.columns
    )
    df = pd.concat((df_in, tot_df), axis=0)

    df.index.rename("experiment", inplace=True)  # give index appropriate name
    return df


@table
def sqrt_datasets_bias_variance_ratio(datasets_bias_variance_ratio):
    """Given `datasets_bias_variance_ratio` take the sqrt and tabulate the
    results. This gives an idea of how
    faithful the uncertainties are in sensible units. As noted in
    `datasets_bias_variance_ratio`, bias/variance is a squared quantity and
    so when considering how much uncertainty has been over or underestimated
    it is more natural to consider sqrt(bias/variance).

    """
    df_in = datasets_bias_variance_ratio
    vals = np.array(df_in.values)  # copy just in case
    vals[:, 1] = np.sqrt(vals[:, 1])
    return pd.DataFrame(
        vals, index=df_in.index, columns=["ndata", "sqrt(bias/variance)"]
    )


@table
def sqrt_experiments_bias_variance_ratio(experiments_bias_variance_ratio):
    """Like sqrt_datasets_bias_variance_ratio except for each experiment.

    """
    return sqrt_datasets_bias_variance_ratio(experiments_bias_variance_ratio)


@table
def total_bias_variance_ratio(
    experiments_bias_variance_ratio, datasets_bias_variance_ratio, experiments
):
    """Combine datasets_bias_variance_ratio and experiments_bias_variance_ratio
    into single table with MultiIndex of experiment and dataset.

    """
    exps_df_in = experiments_bias_variance_ratio.iloc[:-1]  # Handle total separately
    lvs = exps_df_in.index
    # The explicit call to list is because pandas gets confused otherwise
    expanded_index = pd.MultiIndex.from_product((list(lvs), ["Total"]))
    exp_df = exps_df_in.set_index(expanded_index)

    dset_index = pd.MultiIndex.from_arrays(
        [
            [
                str(experiment)
                for experiment in experiments
                for ds in experiment.datasets
            ],
            datasets_bias_variance_ratio.index.values,
        ]
    )
    ds_df = datasets_bias_variance_ratio.set_index(dset_index)
    dfs = []
    for lv in lvs:
        dfs.append(pd.concat((exp_df.loc[lv], ds_df.loc[lv]), copy=False, axis=0))
    total_df = pd.DataFrame(
        experiments_bias_variance_ratio.iloc[[-1]].values,
        columns=exp_df.columns,
        index=["Total"],
    )
    dfs.append(total_df)
    keys = [*lvs, "Total"]
    res = pd.concat(dfs, axis=0, keys=keys)
    return res


@table
def expected_xi_from_bias_variance(sqrt_experiments_bias_variance_ratio):
    """Given the ``sqrt_experiments_bias_variance_ratio`` calculate a predicted
    value of :math:`\\xi_{1 \sigma}` for each experiment. The predicted value is based of
    the assumption that the difference between replica and central prediction
    and the difference between central prediction and underlying prediction are
    both gaussians centered on zero.

    For example, if sqrt(expected bias/expected variance) is 0.5, then we would
    expect xi_{1 sigma} to be given by performing an integral of the
    distribution of

        diffs = (central - underlying predictions)

    over the domain defined by the variance. In this case the sqrt(variance) is
    twice as large as the sqrt(bias) which is the same as integrating a normal
    distribution mean = 0, std = 1 over the interval [-2, 2], given by

        integral = erf(2/sqrt(2))

    where erf is the error function.

    In general the equation is

        integral = erf(sqrt(variance / (2*bias)))

    """
    df_in = sqrt_experiments_bias_variance_ratio
    n_sigma_in_variance = 1 / df_in.values[:, -1, np.newaxis]
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    estimated_integral = special.erf(n_sigma_in_variance / np.sqrt(2))
    return pd.DataFrame(
        np.concatenate((df_in.values[:, 0, np.newaxis], estimated_integral), axis=1),
        index=df_in.index,
        columns=["ndata", r"estimated $\xi_{1\sigma} \,$ from bias/variance"],
    )


@table
def fits_measured_xi(experiments_xi_measured, experiments):
    r"""Tabulate the measure value of \xi_{1\sigma} for each experiment, as
    calculated by experiment_xi. Note that the mean is taken across directions
    of the covariance matrix.

    """
    records = []
    tot_xi = 0
    tot_n = 0
    for exp, xi in zip(experiments, experiments_xi_measured):
        records.append(dict(experiment=str(exp), ndata=len(xi), xi=np.mean(xi)))
        tot_xi += len(xi) * np.mean(xi)
        tot_n += len(xi)
    records.append(dict(experiment="Total", ndata=tot_n, xi=tot_xi / tot_n))
    df = pd.DataFrame.from_records(
        records, index="experiment", columns=("experiment", "ndata", "xi")
    )
    df.columns = ["ndata", r"measured $\xi_{1\sigma}$"]
    return df


@table
def compare_measured_expected_xi(fits_measured_xi, expected_xi_from_bias_variance):
    """Table with measured xi and expected xi from bias/variance for each
    experiment and total. For details on expected xi, see
    expected_xi_from_bias_variance. For more details on measured xi see
    fits_measured_xi.

    """
    # don't want ndata twice
    df = pd.concat(
        (fits_measured_xi, expected_xi_from_bias_variance.iloc[:, 1]), axis=1
    )
    return df


@figure
def plot_dataset_xi(dataset_xi, dataset):
    r"""For a given dataset, plot the value of \xi_{1 \sigma} for each direction
    of the covariance matrix, along with the expected value of \xi_{1 \sigma}
    if the replicas distribution perfectly matches the central distribution
    (0.68). In the legend include the mean across directions.

    """
    fig, ax = plt.subplots()
    ax.plot(
        dataset_xi,
        "*",
        label=r"$\xi_{1\sigma}$ = " + f"{dataset_xi.mean():.2f}, from multifits",
        clip_on=False,
    )
    ax.axhline(
        0.68, linestyle=":", color="k", label=r"$\xi_{1\sigma}$ " + "expected value"
    )
    ax.axhline(
        0.95, linestyle=":", color="r", label=r"$\xi_{2\sigma}$ " + "expected value"
    )
    ax.set_ylim((0, 1))
    ax.set_xlabel("eigenvector index (ascending order)")
    ax.set_title(r"$\xi_{1\sigma}$ for " + str(dataset))
    ax.legend()
    return fig


@figure
def plot_dataset_xi_histogram(dataset_xi, dataset):
    r"""For a given dataset, bin the values of \xi_{1 \sigma} for each direction
    of the covariance matrix, plot as a histogram with a vertical line for the
    expected value: 0.68. In the legend print the mean and standard deviation
    of the distribution.

    """
    fig, ax = plt.subplots()
    ax.hist(
        dataset_xi,
        label=(
            r"$\xi_{1\sigma}$ = "
            + f"{dataset_xi.mean():.2f}, "
            + r"std($\xi_{1\sigma}$) = "
            + f"{dataset_xi.std():.2f}"
        ),
    )
    ax.axvline(
        0.68, linestyle=":", color="k", label=r"$\xi_{1\sigma}$ " + "expected value"
    )
    ax.set_xlim((0, 1))
    ax.set_xlabel(r"$\xi^{i}_{1\sigma}$")
    ax.set_title("Histogram of " + r"$\xi^{i}_{1\sigma}$ for " + str(dataset))
    ax.legend()
    return fig


@figure
def plot_experiment_xi(experiment_xi, experiment):
    """Like plot_dataset_xi except for an experiment.

    """
    return plot_dataset_xi(experiment_xi, experiment)


@figure
def plot_experiment_xi_histogram(experiment_xi, experiment):
    """Like plot_dataset_xi_histogram but for an experiment.

    """
    return plot_dataset_xi_histogram(experiment_xi, experiment)


@table
def dataset_ratio_error_finite_effects(
    bias_variance_resampling_dataset, n_fit_samples, n_replica_samples
):
    """For a single dataset vary number of fits and number of replicas used to perform
    bootstrap sample of expected bias and variance. For each combination of
    n_rep and n_fit tabulate the std deviation across bootstrap samples of

        ratio = bias / variance

    The resulting table gives and approximation of how error varies with
    number of fits and number of replicas for each dataset.

    """
    bias_samples, var_samples = bias_variance_resampling_dataset
    ratio = bias_samples / var_samples
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(ratio.std(axis=2), index=ind, columns=col)


@table
def total_ratio_error_finite_effects(
    bias_variance_resampling_total, n_fit_samples, n_replica_samples
):
    """Like dataset_ratio_relative_error_finite_effects except for the total
    bias / variance (across all data).

    """
    return dataset_ratio_error_finite_effects(
        bias_variance_resampling_total, n_fit_samples, n_replica_samples
    )


@table
def total_ratio_means_finite_effects(
    bias_variance_resampling_total, n_fit_samples, n_replica_samples
):
    """Vary number of fits and number of replicas used to perform
    bootstrap sample of expected bias and variance. For each combination of
    n_rep and n_fit tabulate the the mean across bootstrap samples of

        ratio = total bias / total variance

    which can give context to `total_ratio_relative_error_finite_effects`.

    """
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    bias_total, var_total = bias_variance_resampling_total
    return pd.DataFrame((bias_total / var_total).mean(axis=2), index=ind, columns=col)


@table
def dataset_xi_error_finite_effects(
    xi_resampling_dataset, n_fit_samples, n_replica_samples
):
    """For a single dataset vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the mean of xi across datapoints (note that points
    here refers to points in the basis which diagonalises the covmat) and then
    tabulate the standard deviation of xi_1sigma across bootstrap samples.

    """
    means_xi_table = xi_resampling_dataset.mean(axis=-1)
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(means_xi_table.std(axis=2), index=ind, columns=col)


@table
def dataset_xi_means_finite_effects(
    xi_resampling_dataset, n_fit_samples, n_replica_samples
):
    """For a single dataset vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the mean of xi across datapoints (note that points
    here refers to points in the basis which diagonalises the covmat) and then
    tabulate the mean of xi_1sigma across bootstrap samples.

    """
    means_xi_table = xi_resampling_dataset.mean(axis=-1)
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(means_xi_table.mean(axis=2), index=ind, columns=col)


# NOTE: This action was written when trying to understand the finite size effects
# and is largely redundant.
@table
def dataset_std_xi_error_finite_effects(
    xi_resampling_dataset, n_fit_samples, n_replica_samples
):
    """For a single dataset vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the standard deviation of xi across datapoints
    (note that points here refers to points in the basis which diagonalises the
    covmat) and then tabulate the standard deviation of std(xi_1sigma) across
    bootstrap samples.

    """
    means_xi_table = xi_resampling_dataset.std(axis=-1)
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(means_xi_table.std(axis=2), index=ind, columns=col)


@table
def dataset_std_xi_means_finite_effects(
    xi_resampling_dataset, n_fit_samples, n_replica_samples
):
    """For a single dataset vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the standard deviation of xi across datapoints
    (note that points here refers to points in the basis which diagonalises the
    covmat) and then tabulate the mean of std(xi_1sigma) across
    bootstrap samples.

    """
    means_xi_table = xi_resampling_dataset.std(axis=-1)
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(means_xi_table.mean(axis=2), index=ind, columns=col)


@table
def total_xi_error_finite_effects(total_xi_resample, n_fit_samples, n_replica_samples):
    """For all data vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the mean of xi across datapoints (note that points
    here refers to points in the basis which diagonalises the covmat) and then
    tabulate the standard deviation of xi_1sigma across bootstrap samples.

    """
    return dataset_xi_error_finite_effects(
        total_xi_resample, n_fit_samples, n_replica_samples
    )


@table
def total_xi_means_finite_effects(total_xi_resample, n_fit_samples, n_replica_samples):
    """For all data vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the mean of xi across datapoints (note that points
    here refers to points in the basis which diagonalises the covmat) and then
    tabulate the standard deviation of xi_1sigma across bootstrap samples.

    """
    return dataset_xi_means_finite_effects(
        total_xi_resample, n_fit_samples, n_replica_samples
    )


@table
def total_expected_xi_means_finite_effects(
    total_expected_xi_resample, n_fit_samples, n_replica_samples
):
    """Given the resampled ratio of bias/variance, returns table of mean of
    resampled expected xi across bootstrap samples.

    See `expected_xi_from_bias_variance` for more details on how to calculate
    expected xi.

    """
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(total_expected_xi_resample.mean(axis=2), index=ind, columns=col)


@table
def total_expected_xi_error_finite_effects(
    total_expected_xi_resample, n_fit_samples, n_replica_samples
):
    """Given the resampled ratio of bias/variance, returns table of mean of
    resampled expected xi across bootstrap samples.

    See :py:func:`expected_xi_from_bias_variance` for more details on how to calculate
    expected xi.

    """
    ind = pd.Index(n_replica_samples, name="n_rep samples")
    col = pd.Index(n_fit_samples, name="n_fit samples")
    return pd.DataFrame(total_expected_xi_resample.std(axis=2), index=ind, columns=col)


@table
def total_std_xi_error_finite_effects(
    exps_xi_resample, n_fit_samples, n_replica_samples
):
    """For all data vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the std deviation of xi across datapoints
    (note that points here refers to points in the basis which diagonalises
    the covmat) and then tabulate the mean of std(xi_1sigma) across bootstrap
    samples.

    """
    xi_total = np.concatenate(exps_xi_resample, axis=-1)
    return dataset_std_xi_error_finite_effects(
        xi_total, n_fit_samples, n_replica_samples
    )


@table
def total_std_xi_means_finite_effects(
    exps_xi_resample, n_fit_samples, n_replica_samples
):
    """For all data vary number of fits and number of replicas used to perform
    bootstrap sample of xi. Take the std deviation of xi across datapoints
    (note that points here refers to points in the basis which diagonalises the
    covmat) and then tabulate the standard deviation of std(xi_1sigma) across
    bootstrap samples.

    """
    xi_total = np.concatenate(exps_xi_resample, axis=-1)
    return dataset_std_xi_means_finite_effects(
        xi_total, n_fit_samples, n_replica_samples
    )


@table
def experiments_bootstrap_sqrt_ratio_table(
    experiments_bootstrap_sqrt_ratio, experiments
):
    """Given experiments_bootstrap_sqrt_ratio, which a bootstrap
    resampling of the sqrt(bias/variance) for each experiment and the total
    across all data, tabulate the mean and standard deviation across bootstrap
    samples.

    """
    indices = list(map(str, experiments)) + ["Total"]
    records = []
    for i, exp in enumerate(indices):
        ratio_boot = experiments_bootstrap_sqrt_ratio[i]
        records.append(
            dict(
                experiment=exp,
                mean_ratio=np.mean(ratio_boot, axis=0),
                std_ratio=np.std(ratio_boot, axis=0),
            )
        )
    df = pd.DataFrame.from_records(
        records, index="experiment", columns=("experiment", "mean_ratio", "std_ratio")
    )
    df.columns = [
        "Bootstrap mean sqrt(bias/variance)",
        "Bootstrap std. dev. sqrt(bias/variance)",
    ]
    return df


@table
def experiments_bootstrap_expected_xi_table(
    experiments_bootstrap_expected_xi, experiments
):
    """Tabulate the mean and standard deviation across bootstrap samples of the
    expected xi calculated from the ratio of bias/variance. Returns a table with
    two columns, for the bootstrap mean and standard deviation
    and a row for each experiment plus the total across all experiments.

    """
    df = experiments_bootstrap_sqrt_ratio_table(
        experiments_bootstrap_expected_xi, experiments
    )
    # change the column headers
    df.columns = [
        r"Bootstrap mean expected $\xi_{1\sigma}$ from ratio",
        r"Bootstrap std. dev. expected $\xi_{1\sigma}$ from ratio",
    ]
    return df


@table
def experiments_bootstrap_xi_table(
    experiments_bootstrap_xi, experiments, total_bootstrap_xi
):
    """Tabulate the mean and standard deviation of xi_1sigma across bootstrap
    samples. Note that the mean has already be taken across data points
    (or eigenvector directions in the basis which diagonalises the covariance
    matrix) for each individual bootstrap sample.

    Tabulate the results for each experiment and for the total xi across all data.

    """
    # take mean across data points for each experiment
    xi_1sigma = [np.mean(exp_xi, axis=1) for exp_xi in experiments_bootstrap_xi]
    # take mean across all data
    xi_1sigma.append(np.mean(total_bootstrap_xi, axis=1))
    df = experiments_bootstrap_sqrt_ratio_table(xi_1sigma, experiments)
    df.columns = [
        r"Bootstrap mean $\xi_{1\sigma}$",
        r"Bootstrap std. dev. $\xi_{1\sigma}$",
    ]
    return df


@table
def experiments_bootstrap_xi_comparison(
    experiments_bootstrap_xi_table, experiments_bootstrap_expected_xi_table
):
    """Table comparing the mean and standard deviation across bootstrap samples of
    the measured xi_1sigma and the expected xi_1sigma calculated from
    bias/variance.

    """
    return pd.concat(
        (experiments_bootstrap_xi_table, experiments_bootstrap_expected_xi_table),
        axis=1,
    )
