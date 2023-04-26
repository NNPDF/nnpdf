"""
multiclosure_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for multiclosure estimators in the space of
data.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.special as special
import scipy.stats
import yaml

from reportengine.figure import figure, figuregen
from reportengine.table import table

from validphys.closuretest.multiclosure import expected_dataset_bias_variance
from validphys.loader import Loader

import logging
log = logging.getLogger(__name__)

l = Loader()

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
def plot_data_fits_bias_variance(fits_data_bias_variance, data):
    """Like `plot_dataset_fits_bias_variance` but for all data. Can use
    alongside ``group_dataset_inputs_by_experiment`` to plot for each experiment.

    """
    return plot_dataset_fits_bias_variance(fits_data_bias_variance, data)

@figure
def plot_dataset_fits_sqrt_bias_variance_ratio(fits_dataset_bias_variance, dataset):
    """Given a dataset and a set of closure fits, calculate the sqrt of bias 
    and variance ratio across fits and then plot scatter points.
    
    Parameters
    ----------
    fits_dataset_bias_variance : multiclosure.fits_dataset_bias_variance
                        tuple containing  biases, variances, len(law_th)

    dataset : (DataSetSpec, DataGroupSpec)
        Note that due to the structure of `validphys` this
        function can be overloaded to accept a DataGroupSpec.

    Returns
    -------
    matplotlib fig
    

    """
    biases, variances, _ = fits_dataset_bias_variance
    sqrt_ratio = np.sqrt(biases / variances)

    fig, ax = plt.subplots()
    ax.plot(sqrt_ratio, "*", label = r"$\sqrt{R_{bv}}$ "+f", std. dev. = {np.std(sqrt_ratio):.2f}")
    ax.axhline(np.sqrt(np.mean(biases) / np.mean(variances)), label = f"mean = {np.mean(sqrt_ratio):.2f}")
    ax.set_title(r"$\sqrt{R_{bv}}$ "+f"indicator for {dataset} for each fit")
    ax.set_xlabel("fit index")
    ax.legend()
    return fig

@figure
def plot_data_fits_sqrt_bias_variance_ratio(fits_data_bias_variance, data):
    """
    like `plot_dataset_fits_sqrt_bias_variance_ratio` but for all data.
    """
    return plot_dataset_fits_sqrt_bias_variance_ratio(fits_data_bias_variance,data)
    
@figure
def progressive_sqrt_b_v_ratio_dataset(fits_dataset_bias_variance, dataset):
    """For a set of closure fits, calculate bias and variance across fits on a given dataset.
    Plot the square root ratio between the two quantities as the number of fits increases.
    The progressiv average is calculated as:
    R_{BV} = sqrt(sum_{i = 1}^{n_fits}[bias_i]/sum_{j = 1}^{n_fits}[var_j]))
    with n_fits = [1,2,...,N_tot_fits]

    Parameters
    -------
    fits_dataset_bias_variance : multiclosure.fits_dataset_bias_variance
                        tuple containing  biases, variances, len(law_th)

    dataset : (DataSetSpec, DataGroupSpec)
        Note that due to the structure of `validphys` this
        function can be overloaded to accept a DataGroupSpec.

    Returns
    -------
    matplotlib fig
    

    """
    biases, variances, _ = fits_dataset_bias_variance
    prog_biases = []
    prog_var = []
    for i in range(np.size(biases)):
        prog_biases.append(np.mean(biases[0:i+1]))
        prog_var.append(np.mean(variances[0:i+1]))
    prog_biases = np.asarray(prog_biases)
    prog_var = np.asarray(prog_var)
    fig, ax = plt.subplots()
    ax.plot(np.sqrt(prog_biases/prog_var), "-", label=f"progressive sqrt b/v ratio")

    ax.set_title(f"progressive sqrt b/v ratio for {dataset} for increasing fits")
    ax.set_xlabel("fit index")
    ax.legend()
    return fig

@figure
def plot_total_fits_bias_variance(fits_total_bias_variance):
    """Like `plot_dataset_fits_bias_variance` but for the total bias/variance
    for all data, with the total calculated by summing up contributions from each
    experiment.

    """
    return plot_dataset_fits_bias_variance(fits_total_bias_variance, "all data")

@figure
def progressive_sqrt_b_v_ratio_data(fits_data_bias_variance, data):
    """Like `progressive_sqrt_b_v_ratio` but for all data.

    """
    return progressive_sqrt_b_v_ratio_dataset(fits_data_bias_variance, data)

@figure
def plot_xq2_with_Rbv(dataset_inputs_by_groups_xq2map,datasets_bias_variance_ratio):
    """
    Like validphys.dataplots.plot_xq2 but legend (hue axis) is given by Rbv value on each of
    the datasets
    """
    
    datasets_inp_x2map = dataset_inputs_by_groups_xq2map

    ds_x_q2 = []

    for ds in datasets_inp_x2map:
        ds_x_q2.append({'dataset':str(ds.commondata),'x':ds.fitted[0], 'q2':ds.fitted[1]})
    
    df = pd.DataFrame([(x['dataset'], i, j) for x in ds_x_q2 for i, j in zip(x['x'], x['q2'])], columns=['dataset', 'x', 'q2'])
    rbv_ds = datasets_bias_variance_ratio['sqrt(bias/variance)']
    formatted_rbv_ds = {k: round(v, 2) for k, v in rbv_ds.items()}
    df['rbv'] = df['dataset'].map(dict(formatted_rbv_ds))
    df['hue'] = df['dataset'] + ', ' + df['rbv'].astype(str)

    fig, ax = plt.subplots()
    palette =sns.color_palette("tab10")
    sns.scatterplot(data=df, x="x", y="q2", hue="hue", palette = palette)
    ax.set(xscale="log", yscale="log")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    return fig




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
        records.append(dict(dataset=str(ds), ndata=ndata, ratio=bias / var, sqrt_ratio = np.sqrt(bias/var)))
    df = pd.DataFrame.from_records(
        records, index="dataset", columns=("dataset", "ndata", "ratio","sqrt_ratio")
    )
    df.columns = ["ndata", "bias/variance", "sqrt(bias/variance)"]
    return df

@table
def datasets_bias_variance(datasets_expected_bias_variance, each_dataset):
    """For each dataset calculate the expected bias and expected variance
    across fitsband tabulate the results. Bias and Variance are normalized by number of data points

    Notes
    -----

    This is to check the weight each dataset/process has in the calculation of the complete R_bv ratio. 
    This is because one dataset alone could have a correct B/V=1 but if Bias and Variance are both centered
    around a number != 1 this means that in the calculation of B/V total ratio the specific dataset/
    process is going to have much more weight than the rest

    """
    records = []
    for ds, (bias, var, ndata) in zip(each_dataset, datasets_expected_bias_variance):
        records.append(dict(dataset=str(ds), ndata=ndata, bias=bias, variance=var))
    df = pd.DataFrame.from_records(
        records, index="dataset", columns=("dataset", "ndata", "bias", "variance")
    )
    df.columns = ["ndata", "bias", "variance"]
    return df

@table
def bias_variance_ratio_by_process(datasets_expected_bias_variance, each_dataset):
    """
    Same as `datasets_bias_variance_ratio` but adding another index to table, namely,
    the name of the process
    """
    
    process_type = []

    for ds, (bias, var, ndata) in zip(each_dataset, datasets_expected_bias_variance):
        # open plot file and get process name
        plt_file = ds.commondata.plotfiles[1]
        with open(plt_file,'r') as file:
            card = yaml.safe_load(file)
        prcs = card['nnpdf31_process']

        process_type.append({'process':prcs,'dataset':str(ds),'ndata':ndata,'ratio':bias / var, 
                                        'sqrt_ratio':np.sqrt(bias/var)})

    df = pd.DataFrame(process_type).set_index(['process','dataset'])
    return df

@table
def experiments_bias_variance_ratio(
    experiments_expected_bias_variance, experiments_data, expected_total_bias_variance
):
    """Like datasets_bias_variance_ratio except for each experiment. Also
    calculate and tabulate

        (total expected bias) / (total expected variance)

    where the total refers to summing over all experiments.

    """
    # don't reinvent wheel
    df_in = datasets_bias_variance_ratio(
        experiments_expected_bias_variance, experiments_data
    )

    bias_tot, var_tot, ntotal = expected_total_bias_variance

    tot_df = pd.DataFrame(
        [[ntotal, bias_tot / var_tot, np.sqrt(bias_tot / var_tot)]], index=["Total"], columns=df_in.columns
    )
    df = pd.concat((df_in, tot_df), axis=0)
    
    df.index.rename("experiment", inplace=True)  # give index appropriate name
    return df

@table
def experiments_bias_variance(
    experiments_expected_bias_variance, experiments_data, expected_total_bias_variance
):
    """Like datasets_bias_variance_ratio except for each experiment. Also
    calculate and tabulate

        (total expected bias) / (total expected variance)

    where the total refers to summing over all experiments.

    """
    # don't reinvent wheel
    df_in = datasets_bias_variance(
        experiments_expected_bias_variance, experiments_data
    )

    bias_tot, var_tot, ntotal = expected_total_bias_variance

    tot_df = pd.DataFrame(
        [[ntotal, bias_tot/ntotal, var_tot/ntotal]], index=["Total"], columns=df_in.columns
    )
    df = pd.concat((df_in, tot_df), axis=0)
    
    df.index.rename("experiment", inplace=True)  # give index appropriate name
    return df

@table
def experiments_bias_variance_table(
    experiments_expected_bias_variance,
    group_dataset_inputs_by_experiment,
    expected_total_bias_variance,
):
    """Tabulate the values of bias and variance for each experiment as well
    as the sqrt ratio of the two as in
    :py:func`sqrt_experiments_bias_variance_ratio`. Used as a performance
    indicator.

    """
    records = []
    for exp, (bias, var, ndata) in zip(
        group_dataset_inputs_by_experiment,
        experiments_expected_bias_variance
    ):
        records.append(dict(
            experiment=exp["group_name"],
            ndata=ndata,
            bias=bias/ndata,
            variance=var/ndata,
            sqrt_ratio=np.sqrt(bias/var)
        ))

    bias_tot, var_tot, ntotal = expected_total_bias_variance

    records.append(dict(
        experiment="Total",
        ndata=ntotal,
        bias=bias_tot/ntotal,
        variance=var_tot/ntotal,
        sqrt_ratio=np.sqrt(bias_tot/var_tot)
    ))
    df = pd.DataFrame.from_records(records, index="experiment")
    df.columns = [
        "ndata",
        "bias",
        "variance",
        "sqrt(bias/variance)"
    ]
    return df


@table
def total_bias_variance_ratio(
    experiments_bias_variance_ratio, datasets_bias_variance_ratio, experiments_data
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
                for experiment in experiments_data
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
def fits_measured_xi(experiments_xi_measured, experiments_data):
    r"""Tabulate the measure value of \xi_{1\sigma} for each experiment, as
    calculated by data_xi (collected over experiments). Note that the mean is
    taken across eigenvectors of the covariance matrix.

    """
    records = []
    tot_xi = 0
    tot_n = 0
    for exp, xi in zip(experiments_data, experiments_xi_measured):
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
    r"""For a given dataset, plot the value of \xi_{1 \sigma} for each eigenvector
    of the covariance matrix, along with the expected value of \xi_{1 \sigma}
    if the replicas distribution perfectly matches the central distribution
    (0.68). In the legend include the mean across eigenvectors.

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
    r"""For a given dataset, bin the values of \xi_{1 \sigma} for each eigenvector
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
def plot_data_xi(data_xi, data):
    """Like plot_dataset_xi except for all data. Can be used alongside
    ``group_dataset_inputs_by_experiment`` to plot for each experiment.

    """
    return plot_dataset_xi(data_xi, data)


@figure
def plot_data_xi_histogram(data_xi, data):
    """Like plot_dataset_xi_histogram but for all data. Can be used alongside
    ``group_dataset_inputs_by_experiment`` to plot for each experiment.

    """
    return plot_dataset_xi_histogram(data_xi, data)


@figure
def plot_data_central_diff_histogram(experiments_replica_central_diff):
    """Histogram of the difference between central prediction
    and underlying law prediction normalised by the corresponding replica
    standard deviation by concatenating the difference across all data. plot a
    scaled gaussian for reference. Total xi is the number of central differences
    which fall within the 1-sigma confidence interval of the scaled gaussian.

    """
    scaled_diffs = np.concatenate([
        (central_diff / sigma).flatten()
        for sigma, central_diff
        in experiments_replica_central_diff
    ])
    fig, ax = plt.subplots()
    ax.hist(
        scaled_diffs, bins=50, density=True, label="Central prediction distribution"
    )
    xlim = (-5, 5)
    ax.set_xlim(xlim)

    x = np.linspace(*xlim, 100)
    ax.plot(
        x,
        scipy.stats.norm.pdf(x),
        "-k",
        label="Normal distribution",
    )
    ax.legend()
    ax.set_xlabel("Difference to underlying prediction")
    return fig



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
    experiments_bootstrap_sqrt_ratio, experiments_data
):
    """Given experiments_bootstrap_sqrt_ratio, which a bootstrap
    resampling of the sqrt(bias/variance) for each experiment and the total
    across all data, tabulate the mean and standard deviation across bootstrap
    samples.

    """
    indices = list(map(str, experiments_data)) + ["Total"]
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
def groups_bootstrap_sqrt_ratio_table(
    groups_bootstrap_sqrt_ratio, groups_data
):
    """Like :py:func:`experiments_bootstrap_sqrt_ratio_table` but for
    metadata groups.
    """
    df = experiments_bootstrap_sqrt_ratio_table(
        groups_bootstrap_sqrt_ratio, groups_data
    )
    idx = df.index.rename("group")
    return df.set_index(idx)


@table
def experiments_bootstrap_expected_xi_table(
    experiments_bootstrap_expected_xi, experiments_data
):
    """Tabulate the mean and standard deviation across bootstrap samples of the
    expected xi calculated from the ratio of bias/variance. Returns a table with
    two columns, for the bootstrap mean and standard deviation
    and a row for each experiment plus the total across all experiments.

    """
    df = experiments_bootstrap_sqrt_ratio_table(
        experiments_bootstrap_expected_xi, experiments_data
    )
    # change the column headers
    df.columns = [
        r"Bootstrap mean expected $\xi_{1\sigma}$ from ratio",
        r"Bootstrap std. dev. expected $\xi_{1\sigma}$ from ratio",
    ]
    return df


@table
def groups_bootstrap_expected_xi_table(groups_bootstrap_expected_xi, groups_data):
    """Like :py:func:`experiments_bootstrap_expected_xi_table` but for metadata
    groups.
    """
    df = experiments_bootstrap_expected_xi_table(
        groups_bootstrap_expected_xi, groups_data)
    idx = df.index.rename("group")
    return df.set_index(idx)

@table
def experiments_bootstrap_xi_table(
    experiments_bootstrap_xi, experiments_data, total_bootstrap_xi
):
    """Tabulate the mean and standard deviation of xi_1sigma across bootstrap
    samples. Note that the mean has already be taken across data points
    (or eigenvectors in the basis which diagonalises the covariance
    matrix) for each individual bootstrap sample.

    Tabulate the results for each experiment and for the total xi across all data.

    """
    # take mean across data points for each experiment
    xi_1sigma = [np.mean(exp_xi, axis=1) for exp_xi in experiments_bootstrap_xi]
    # take mean across all data
    xi_1sigma.append(np.mean(total_bootstrap_xi, axis=1))
    df = experiments_bootstrap_sqrt_ratio_table(xi_1sigma, experiments_data)
    df.columns = [
        r"Bootstrap mean $\xi_{1\sigma}$",
        r"Bootstrap std. dev. $\xi_{1\sigma}$",
    ]
    return df


@table
def groups_bootstrap_xi_table(
    groups_bootstrap_xi, groups_data, total_bootstrap_xi
):
    """Like :py:func:`experiments_bootstrap_xi_table` but for metadata groups."""
    df = experiments_bootstrap_xi_table(
        groups_bootstrap_xi, groups_data, total_bootstrap_xi)
    idx = df.index.rename("group")
    return df.set_index(idx)


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


@table
def groups_bootstrap_xi_comparison(
    groups_bootstrap_xi_table, groups_bootstrap_expected_xi_table
):
    """Like :py:func:`experiments_bootstrap_xi_comparison` but for metadata
    groups.
    """
    return experiments_bootstrap_xi_comparison(
        groups_bootstrap_xi_table, groups_bootstrap_expected_xi_table
    )


@figuregen
def plot_experiments_sqrt_ratio_bootstrap_distribution(
    experiments_bootstrap_sqrt_ratio, experiments_data
):
    """Plots a histogram for each experiment and the total, showing the
    distribution of bootstrap samples. Takes the mean and std deviation of the
    bootstrap sample and plots the corresponding scaled normal distribution
    for comparison. The limits are set to be +/- 3 std deviations of the mean.

    """
    # experiments_bootstrap_sqrt_ratio includes total. str(exp) is only used to
    # generate title, so appending string is fine.
    for sqrt_ratio_sample, exp in zip(
        experiments_bootstrap_sqrt_ratio, experiments_data + ["Total"]
    ):
        fig, ax = plt.subplots()
        ax.hist(sqrt_ratio_sample, bins=20, density=True)
        mean = np.mean(sqrt_ratio_sample)
        std = np.std(sqrt_ratio_sample)

        xlim = (mean - 3 * std, mean + 3 * std)
        ax.set_xlim(xlim)

        x = np.linspace(*xlim, 100)
        ax.plot(
            x,
            scipy.stats.norm.pdf(x, mean, std),
            "-r",
            label=f"Corresponding normal distribution: mean = {mean:.2g}, std = {std:.2g}",
        )
        ax.legend()
        ax.set_title(f"Bootstrap distribution of sqrt(bias/variance) for {exp}")
        ax.set_xlabel("Sqrt(bias/variance)")
        yield fig


@figuregen
def plot_experiments_xi_bootstrap_distribution(
    experiments_bootstrap_xi, total_bootstrap_xi, experiments_data
):
    """Similar to :py:func:`plot_sqrt_ratio_bootstrap_distribution` except plots
    the bootstrap distribution of xi_1sigma, along with a corresponding
    scaled gaussian, for each experiment (and for all data).

    """
    # take mean across data points for each experiment
    xi_1sigma = [np.mean(exp_xi, axis=1) for exp_xi in experiments_bootstrap_xi]
    # take mean across all data
    xi_1sigma.append(np.mean(total_bootstrap_xi, axis=1))
    # use plotting function from above
    xi_plots = plot_experiments_sqrt_ratio_bootstrap_distribution(
        xi_1sigma, experiments_data
    )
    # Update the title and x label on each plot to reflect that we're plotting
    # \xi_1sigma, don't forget Total plot.
    for fig, exp in zip(xi_plots, experiments_data + ["Total"]):
        ax = fig.gca()
        ax.set_title(r"Bootstrap distribution of $\xi_{1\sigma}$ for " + str(exp))
        ax.set_xlabel(r"$\xi_{1\sigma}$")
        yield fig

@figuregen
def plot_bias_variance_distributions(
    experiments_fits_bias_replicas_variance_samples,
    group_dataset_inputs_by_experiment
):
    """For each experiment, plot the distribution across fits of bias
    and the distribution across fits and replicas of

    fit_rep_var = (E[g] - g)_i inv(cov)_ij (E[g] - g)_j

    where g is the replica prediction for fit l, replica k and E[g] is the
    mean across replicas of g for fit l.

    """
    for (exp_biases, exp_vars, _), group_spec in zip(
            experiments_fits_bias_replicas_variance_samples,
            group_dataset_inputs_by_experiment
        ):
        fig, ax = plt.subplots()
        labels = [
            "fits bias distribution",
            "replicas variance distribution",
        ]
        ax.hist(
            [exp_biases, exp_vars],
            density=True,
            label=labels
        )
        ax.legend()
        ax.set_title(
            f"Bias and variance distributions for {group_spec['group_name']}."
        )
        yield fig
    total_bias, total_var, _ = np.sum(
        experiments_fits_bias_replicas_variance_samples,
        axis=0
    )
    fig, ax = plt.subplots()
    ax.hist(
        [total_bias, total_var],
        density=True,
        label=labels
    )
    ax.legend()
    ax.set_title("Total bias and variance distributions.")
    yield fig
