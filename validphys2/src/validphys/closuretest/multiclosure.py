"""
closuretest/multiclosure.py

containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit
"""

import numpy as np
import scipy.linalg as la
import scipy.special as special
import pandas as pd
import matplotlib.pyplot as plt

from reportengine import collect
from reportengine.table import table
from reportengine.figure import figure

from validphys.results import ThPredictionsResult
from validphys.calcutils import calc_chi2
from validphys.core import DataSetSpec, Experiment

def internal_multiclosure_dataset_loader(
    dataset, fits_pdf, multiclosure_underlyinglaw):
    """internal function for loading multiple theory predictions for a given
    experiment, a single covariance matrix using underlying law as t0pdf for use
    with multiclosure statistical estimators. Avoiding memory issues from caching
    experiment load

    Returns
    -------
    multiclosure_results: tuple
        a tuple of length 4 containing all necessary dependencies of multiclosure
        statistical estimators in order:

            closure fits theory predictions,
            underlying law theory predictions,
            covariance matrix, 
            sqrt covariance matrix

    TODO: deprecate this at some point
    """
    if isinstance(dataset, DataSetSpec):
        data = dataset.load() #just use internal loader
    else:
        # it's an experiment, so copy loading of ExperimentSpec.core without
        # caching
        sets = []
        for ds in dataset.datasets:
            loaded_data = ds.load()
            sets.append(loaded_data)
        data = Experiment(sets, dataset.name)

    fits_dataset_predictions = [
        ThPredictionsResult.from_convolution(pdf, dataset, loaded_data=data)
        for pdf in fits_pdf
    ]
    fits_underlying_predictions = ThPredictionsResult.from_convolution(
        multiclosure_underlyinglaw, dataset, loaded_data=data)

    # copy data to make t0 cov
    loaded_data = type(data)(data)
    loaded_data.SetT0(multiclosure_underlyinglaw.load_t0())
    covmat = loaded_data.get_covmat()
    sqrt_covmat = la.cholesky(covmat, lower=True)
    # TODO: support covmat reg and theory covariance matrix
    # possibly make this a named tuple
    return (
        fits_dataset_predictions,
        fits_underlying_predictions,
        covmat,
        sqrt_covmat
    )

def internal_multiclosure_experiment_loader(
    experiment, fits_pdf, multiclosure_underlyinglaw
):
    """Like `internal_multiclosure_dataset_loader` except for an experiment"""
    return internal_multiclosure_dataset_loader(
        experiment, fits_pdf, multiclosure_underlyinglaw)

def expected_bias_dataset(
    internal_multiclosure_dataset_loader):
    """For a set of closure fits, calculate the mean bias across fits for
    a dataset, bias is the chi2 between central prediction and underlying law.
    For more information on bias, see closuretest.bias_dataset.

    The fits should each have the same underlying law and t0pdf, but have
    different filterseeds, so that the level 1 shift is different.

    """
    closures_th, law_th, _, sqrtcov = internal_multiclosure_dataset_loader
    centrals = np.asarray(
        [th.central_value for th in closures_th]
    )
    # place bins on first axis
    diffs = (
        law_th.central_value[:, np.newaxis] -
        centrals.T
    )
    biases = calc_chi2(
        sqrtcov, diffs)
    return np.mean(biases), len(law_th)

def expected_variance_dataset(internal_multiclosure_dataset_loader):
    """Given multiple closure fits, calculate the mean variance across fits
    for a given dataset, variance is the spread of replicas in units of the
    covariance, for more information see closuretest.variance_dataset

    As with expected_bias_dataset, the fits should each have the same underlying
    law and t0pdf, but have different filterseeds.

    """
    closures_th, law_th, _, sqrtcov = internal_multiclosure_dataset_loader
    # object is fits, bins, replicas
    reps = np.asarray([th._rawdata for th in closures_th])
    variances = []
    # this seems slow but breaks for datasets with single data point otherwise
    for i in range(reps.shape[0]):
        diffs = reps[i, :, :] - reps[i, :, :].mean(axis=1, keepdims=True)
        variances.append(np.mean(calc_chi2(sqrtcov, diffs)))
    # assume all data same length
    return np.mean(variances), len(law_th)

def expected_bias_experiment(internal_multiclosure_experiment_loader):
    """Like expected_bias_dataset but for whole experiment"""
    return expected_bias_dataset(internal_multiclosure_experiment_loader)

def expected_variance_experiment(internal_multiclosure_experiment_loader):
    """Like expected_variance_dataset but for whole experiment"""
    return expected_variance_dataset(internal_multiclosure_experiment_loader)

datasets_expected_bias = collect(
    "expected_bias_dataset", ("experiments", "experiment")
)
datasets_expected_variance = collect(
    "expected_variance_dataset", ("experiments", "experiment")
)

#TODO: check that this hasn't been implemented somewhere else at point of merge
experiments_datasets = collect("dataset", ("experiments", "experiment"))

@table
def datasets_bias_variance_ratio(
    datasets_expected_bias, datasets_expected_variance, experiments_datasets
):
    """For each dataset calculate the expected bias and expected variance
    across fits

        ratio = expected bias / expected variance

    and tabulate the results.

    This gives an idea of how faithful uncertainties are for a set of
    datasets.

    """
    records = []
    for ds, (bias, ndata), (var, _) in zip(
        experiments_datasets, datasets_expected_bias, datasets_expected_variance
    ):
        records.append(dict(
            dataset=str(ds),
            ndata=ndata,
            ratio=bias/var
        ))
    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "ndata", "ratio")
    )
    df.columns = ["ndata", "bias/variance"]
    return df

experiments_expected_bias = collect(
    "expected_bias_experiment", ("experiments",)
)
experiments_expected_variance = collect(
    "expected_variance_experiment", ("experiments",)
)

@table
def experiments_bias_variance_ratio(
    experiments_expected_bias, experiments_expected_variance, experiments):
    """Like datasets_bias_variance_ratio except for each experiment. Also
    calculate and tabulate

        total expected bias / total expected variance.

    """
    # don't reinvent wheel
    df_in = datasets_bias_variance_ratio(
        experiments_expected_bias, experiments_expected_variance, experiments)

    bias_tot = np.sum([bias for (bias, _) in experiments_expected_bias])
    var_tot = np.sum([var for (var, _) in experiments_expected_variance])
    ntotal = np.sum(df_in['ndata'].values)

    tot_df = pd.DataFrame(
        [[ntotal, bias_tot/var_tot]], index=["Total"], columns=df_in.columns)
    df = pd.concat((df_in, tot_df), axis=0)

    df.index.rename("experiment", inplace=True) # give index appropriate name
    return df

@table
def sqrt_datasets_bias_variance_ratio(datasets_bias_variance_ratio):
    """For each of the tabulated bias/variance from
    `datasets_bias_variance_ratio` take the sqrt, which gives an idea of how
    faithful the uncertainties are in units of the standard deviation.

    """
    df_in = datasets_bias_variance_ratio
    vals = np.array(df_in.values) # copy just in case
    vals[:, 1] = np.sqrt(vals[:, 1])
    return pd.DataFrame(
        vals,
        index=df_in.index,
        columns=["ndata", "sqrt(bias/variance)"]
    )

@table
def sqrt_experiments_bias_variance_ratio(experiments_bias_variance_ratio):
    """Like sqrt_datasets_bias_variance_ratio except for each experiment
    """
    return sqrt_datasets_bias_variance_ratio(experiments_bias_variance_ratio)

@table
def total_bias_variance_ratio(
    experiments_bias_variance_ratio, datasets_bias_variance_ratio, experiments):
    """Combine datasets_bias_variance_ratio and experiments_bias_variance_ratio
    into single table with multiindex of experiment and dataset
    """
    exps_df_in = experiments_bias_variance_ratio.iloc[:-1] # Handle total seperate
    lvs = exps_df_in.index
    #The explicit call to list is because pandas gets confused otherwise
    expanded_index = pd.MultiIndex.from_product(
        (list(lvs), ["Total"]),
    )
    exp_df = exps_df_in.set_index(expanded_index)

    dset_index = pd.MultiIndex.from_arrays(
        [
            [str(experiment) for experiment in experiments for ds in experiment.datasets],
            datasets_bias_variance_ratio.index.values
        ],
    )
    ds_df = datasets_bias_variance_ratio.set_index(dset_index)
    dfs = []
    for lv in lvs:
        dfs.append(
            pd.concat((exp_df.loc[lv], ds_df.loc[lv]), copy=False, axis=0))
    total_df = pd.DataFrame(
        experiments_bias_variance_ratio.iloc[[-1]].values,
        columns=exp_df.columns,
        index=['Total'],
    )
    dfs.append(total_df)
    keys = [*lvs, 'Total']
    res = pd.concat(dfs, axis=0, keys=keys)
    return res

@table
def expected_xi_from_bias_variance(sqrt_experiments_bias_variance_ratio):
    """Given the `sqrt_experiments_bias_variance_ratio` calculate a predicted
    value of xi_{1 sigma} for each experiment. The predicted value is based of
    the assumption that the difference between replica and central prediction
    and the difference between central prediction and underlying prediction are
    both gaussians centered on zero.

    For example, if sqrt(expected bias/expected variance) is 0.5, then we would
    expect xi_{1 sigma} to be given by performing integral of the distribution
    of

        diffs = (central - underlying predictions)

    over the domain defined by the variance. In this case the sqrt(variance) is
    twice as large as the sqrt(bias) which is the same as integrating a normal
    distribution mean = 0, std = 1 over the interval [-2, 2], given by

        integral = erf(2/sqrt(2))

    where erf is the error function.

    In general the equations is

        integral = erf(sqrt(variance / (2*bias)))

    """
    df_in = sqrt_experiments_bias_variance_ratio
    n_sigma_in_variance = 1/df_in.values[:, -1, np.newaxis]
    # pylint can't find erf here, disable error in this function
    #pylint: disable=no-member
    estimated_integral = special.erf(n_sigma_in_variance/np.sqrt(2))
    return pd.DataFrame(
        np.concatenate((df_in.values[:, 0, np.newaxis], estimated_integral), axis=1),
        index=df_in.index,
        columns=["ndata", r"$\xi \,$ from ratio"]
    )

def dataset_xi(
    internal_multiclosure_dataset_loader):
    r"""For a given dataset calculate sigma, the RMS difference between
    replica predictions and central predictions, and delta, the difference
    between the central prediction and the underlying prediction.

    The differences are calculated in the basis which would diagonalise the
    dataset's covariance matrix.

    Then the Indictor function is evaluated for elementwise for sigma and delta

        I_{[-\sigma_j, \sigma_j]}(\delta_j)

    which is 1 when |\delta_j| < \sigma_j and 0 otherwise. Finally, take the
    mean across fits

    """
    closures_th, law_th, covmat, _ = internal_multiclosure_dataset_loader
    replicas = np.asarray([th._rawdata for th in closures_th])
    centrals = np.mean(replicas, axis=-1)
    underlying = law_th.central_value

    _, e_vec = la.eigh(covmat)

    central_diff = centrals - underlying[np.newaxis, :]
    var_diff_sqrt = (centrals[:, :, np.newaxis] - replicas)

    # project into basis which diagonalises covariance matrix
    var_diff_sqrt = e_vec.T @ var_diff_sqrt.transpose(2, 1, 0)
    central_diff = e_vec.T @ central_diff.T

    var_diff = (var_diff_sqrt)**2
    sigma = np.sqrt(var_diff.mean(axis=0)) # sigma is always positive
    in_1_sigma = np.array(abs(central_diff) < sigma, dtype=int)
    # mean across fits
    return in_1_sigma.mean(axis=1)

def experiment_xi(internal_multiclosure_experiment_loader):
    """Like dataset_xi but for whole experiment"""
    return dataset_xi(internal_multiclosure_experiment_loader)

experiments_xi_measured = collect("experiment_xi", ("experiments",))

@table
def fits_measured_xi(experiments_xi_measured, experiments):
    r"""Tabulate the measure value of \xi_{1\sigma} for each experiment, as
    calculated by experiment_xi, note that the mean is taken across directions
    of the covariance matrix

    """
    records = []
    for exp, xi in zip(
        experiments, experiments_xi_measured
    ):
        records.append(dict(
            experiment=str(exp),
            ndata=len(xi),
            xi=np.mean(xi)
        ))
    df = pd.DataFrame.from_records(
        records,
        index="experiment",
        columns=("experiment", "ndata", "xi")
    )
    df.columns = ["ndata", r"measured $\xi_{1\sigma}$"]
    return df

@figure
def plot_dataset_xi(dataset_xi, dataset):
    r"""For a given dataset, plot the value of \xi_{1 \sigma} for each direction
    of the covariance matrix, along with the expected value of \xi_{1 \sigma}
    if the replicas distribution perfectly matches the central disribution
    (0.68). In the legend include the mean across directions

    """
    fig, ax = plt.subplots()
    ax.plot(
        dataset_xi,
        "*",
        label=r"$\xi_{1\sigma}$ = "+f"{dataset_xi.mean():.2f}, from multifits"
    )
    ax.axhline(0.68, linestyle=":", color="k", label=r"$\xi_{1\sigma}$ "+"expected value")
    ax.axhline(0.95, linestyle=":", color="r", label=r"$\xi_{2\sigma}$"+"expected value")
    ax.set_ylim((0, 1))
    ax.set_xlabel("eigenvector index (ascending order)")
    ax.set_title(r"$\xi_{1\sigma}$ for "+str(dataset))
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
            r"$\xi_{1\sigma}$ = " +
            f"{dataset_xi.mean():.2f}, " +
            r"std($\xi_{1\sigma}$) = " +
            f"{dataset_xi.std():.2f}"
        )
    )
    ax.axvline(0.68, linestyle=":", color="k", label=r"$\xi_{1\sigma}$ "+"expected value")
    ax.set_xlim((0, 1))
    ax.set_xlabel(r"$\xi^{i}_{1\sigma}$")
    ax.set_title("Histogram of " + r"$\xi^{i}_{1\sigma}$ for " + str(dataset))
    ax.legend()
    return fig

@figure
def plot_experiment_xi(experiment_xi, experiment):
    """Like plot_dataset_xi except for an experiment"""
    return plot_dataset_xi(experiment_xi, experiment)

@figure
def plot_experiment_xi_histogram(experiment_xi, experiment):
    """Like plot_dataset_xi_histogram but for an experiment"""
    return plot_dataset_xi_histogram(experiment_xi, experiment)
