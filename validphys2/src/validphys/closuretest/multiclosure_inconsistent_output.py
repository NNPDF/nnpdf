"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

TODO: 
- move some of the functions to multiclosure output
- move the module in closuretest
"""

import pandas as pd
import numpy as np

from reportengine.figure import figuregen, figure
from reportengine.table import table
from reportengine import collect
from scipy.stats import norm

from validphys import plotutils
from validphys.closuretest.multiclosure import regularized_multiclosure_dataset_loader


@table
def table_bias_datasets(bias_datasets, each_dataset):
    """
    TODO: should stay in multiclosure_output


    Compute the bias and sqrt bias and associated errors for each dataset
    and return a DataFrame with the results.

    Parameters
    ----------

    bias_datasets: list
        List of tuples containing the values of bias for each dataset.

    each_dataset: list
        List of validphys.core.DataSetSpec

    Returns
    -------
    pd.DataFrame
        DataFrame containing the bias, variance, ratio and sqrt(ratio) for each dataset
    """
    records = []
    for (biases, n_comp), ds in zip(bias_datasets, each_dataset):
        bias = np.mean(biases)
        std_bias = np.std(biases)

        sqrt_bias = np.sqrt(bias)
        # use gaussian uncertainty propagation
        err_sqrt_bias = 0.5 * std_bias / sqrt_bias

        records.append(
            {
                "dataset": str(ds),
                "dof": n_comp,
                "bias": bias,
                "err bias": std_bias,
                "sqrt bias": sqrt_bias,
                "err sqrt bias": err_sqrt_bias,
            }
        )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "dof", "bias", "err bias", "sqrt bias", "err sqrt bias"),
    )

    return df


@table
def table_bias_data(bias_data):
    """
    TODO: should stay in multiclosure_output

    Same as table_bias_datasets but for all the data, meaning that
    the correlations between the datasets are taken into account.

    Parameters
    ----------
    bias_data: list
        Same of bias_dataset but for all the data

    Returns
    -------
    pd.DataFrame
        DataFrame containing the bias, variance, ratio and sqrt(ratio) for each dataset
    """
    records = []

    biases_tot, n_comp_tot = bias_data
    bias_tot = np.mean(biases_tot)
    std_bias_tot = np.std(biases_tot)

    sqrt_bias_tot = np.sqrt(bias_tot)
    # use gaussian uncertainty propagation
    err_sqrt_bias_tot = 0.5 * std_bias_tot / sqrt_bias_tot

    records.append(
        {
            "dataset": "Total",
            "dof": n_comp_tot,
            "bias": bias_tot,
            "err bias": std_bias_tot,
            "sqrt bias": sqrt_bias_tot,
            "err sqrt bias": err_sqrt_bias_tot,
        }
    )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "dof", "bias", "err bias", "sqrt bias", "err sqrt bias"),
    )
    return df


lambdavalues_table_bias_variance_datasets = collect(
    "table_bias_variance_datasets", ("lambdavalues",)
)


@table
def bootstrapped_table_bias_datasets(bootstrapped_bias_datasets):
    """
    TODO: should stay in multiclosure_output


    Compute the bias, variance, ratio and sqrt(ratio) for each dataset
    and return a DataFrame with the results.
    Uncertainty on ratio and sqrt ratio is computed by Gaussian error propagation
    of the bootstrap uncertainty on bias and variance.
    """
    records = []
    for boot_ds in bootstrapped_bias_datasets:
        df = boot_ds
        mean_bias = df["bias"].mean()

        # gaussian error propagation for the ratio of the means uncertainty
        # only consider bias as source of uncertainty for the ratio (variance is almost constant)
        bootstrap_unc = np.std(df["bias"])
        sqrt_bias = np.mean(np.sqrt(df["bias"]))

        # gaussian error propagation for the sqrt of the ratio
        bootstrap_unc_sqrt = np.std(np.sqrt(df["bias"]))

        records.append(
            {
                "dataset": df["dataset"].iloc[0],
                "mean_dof": df.n_comp.mean(),
                "bias": mean_bias,
                "err bias": bootstrap_unc,
                "sqrt bias": sqrt_bias,
                "err sqrt bias": bootstrap_unc_sqrt,
            }
        )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "mean_dof", "bias", "err bias", "sqrt bias", "err sqrt bias"),
    )
    return df


"""
TODO
"""
lambdavalues_bootstrapped_table_bias_datasets = collect(
    "bootstrapped_table_bias_datasets", ("lambdavalues",)
)


@table
def bootstrapped_table_bias_data(bootstrapped_bias_data):
    """
    TODO: should be in multiclosure_output


    Compute the bias, sqrt bias and their bootstrap errors for a DataGroup
    and return a DataFrame with the results.
    """
    records = []

    df = bootstrapped_bias_data

    bias = df["bias"].mean()
    sqrt_bias = np.mean(np.sqrt(df["bias"].values))

    boot_err_bias = np.std(df["bias"])
    boot_err_sqrt_bias = np.std(np.sqrt(df["bias"]))

    records.append(
        {
            "dataset": "Total",
            "mean_dof": df.n_comp.mean(),
            "bias": bias,
            "err bias": boot_err_bias,
            "sqrt bias": sqrt_bias,
            "err sqrt bias": boot_err_sqrt_bias,
        }
    )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "mean_dof", "bias", "err bias", "sqrt bias", "err sqrt bias"),
    )
    return df


"""
TODO
"""
lambdavalues_bootstrapped_table_bias_data = collect(
    "bootstrapped_table_bias_data", ("lambdavalues",)
)


@figuregen
def plot_lambdavalues_bias_values(
    lambdavalues_bootstrapped_table_bias_datasets, lambdavalues, each_dataset
):
    """
    Plot sqrt of bias and its bootstrap uncertainty as a function of lambda for each dataset.

    Parameters
    ----------
    lambdavalues_bootstrapped_table_bias_datasets: list
        list of data frames computed as per table_bias_datasets.

    lambdavalues: list
        list specified in multiclosure_analysis.yaml

    each_dataset: list
        list of datasets

    Yields
    ------
    figure
    """

    for ds in each_dataset:
        fig, ax = plotutils.subplots()
        for i, lambdavalue in enumerate(lambdavalues):
            df = lambdavalues_bootstrapped_table_bias_datasets[i]
            df = df[df.index == str(ds)]

            ax.errorbar(
                lambdavalue["lambda_value"],
                df["sqrt bias"].values,
                yerr=df["err sqrt bias"].values,
                color="blue",
                fmt='o',
            )
            ax.hlines(1, xmin=0, xmax=1.0, color="red", linestyle="--")
            ax.set_ylabel(r"$R_{b}$")
            ax.set_xlabel(r"$\lambda$")

        ax.set_title(f"{ds.commondata.metadata.plotting.dataset_label}")

        yield fig


@figure
def plot_lambdavalues_bias_values_full_data(
    lambdavalues_bootstrapped_table_bias_data, lambdavalues
):
    """
    Plot sqrt of bias and its bootstrap uncertainty as a function of lambda for the full dataset.

    Parameters
    ----------
    lambdavalues_bootstrapped_table_bias_data: list
        list of data frames computed as per table_bias_data.

    lambdavalues: list
        list specified in multiclosure_analysis.yaml

    Returns
    -------
    figure
    """
    fig, ax = plotutils.subplots()
    for i, lambdavalue in enumerate(lambdavalues):
        df = lambdavalues_bootstrapped_table_bias_data[i]

        ax.errorbar(
            lambdavalue["lambda_value"],
            df["sqrt bias"].values,
            yerr=df["err sqrt bias"].values,
            color="blue",
            fmt='o',
        )
        ax.hlines(1, xmin=0, xmax=1.0, color="red", linestyle="--")
        ax.set_ylabel(r"$R_{b}$")
        ax.set_xlabel(r"$\lambda$")

    ax.set_title("Full data")

    return fig


internal_multiclosure_data_collected_loader = collect("multiclosure_dataset_loader", ("data",))


@figuregen
def plot_l2_condition_number(
    each_dataset, internal_multiclosure_data_collected_loader, evr_min=0.90, evr_max=0.995, evr_n=20
):
    """
    Plot the L2 condition number of the covariance matrix as a function of the explained variance ratio.
    The plot gives an idea of the stability of the covariance matrix as a function of the
    exaplained variance ratio and hence the number of principal components used to reduce the dimensionality.

    The ideal explained variance ratio is chosen based on a threshold L2 condition number, in general this
    threshold number (and the derived explained variance ratio) should be chosen so that

    relative error in output (inverse covmat) <= relative error in input (covmat) * condition number
    Note that in a closure test the relative error in the covariance matrix is very small and only numerical.

    Parameters
    ----------
    each_dataset : list
        List of datasets

    multiclosure_data_loader: list
        list of multiclosure_dataset_loader objects


    Yields
    ------
    fig
        Figure object
    """

    # Explained variance ratio range
    evr_range = np.linspace(evr_min, evr_max, evr_n)

    for multiclosure_dataset_loader, ds in zip(
        internal_multiclosure_data_collected_loader, each_dataset
    ):
        l2_cond = []
        dof = []

        for evr in evr_range:

            pca_loader = regularized_multiclosure_dataset_loader(
                multiclosure_dataset_loader, explained_variance_ratio=evr
            )

            # if the number of principal components is 1 then the covariance matrix is a scalar
            # and the condition number is not computed
            if pca_loader.n_comp == 1:
                l2_cond.append(np.nan)
            else:
                l2_cond.append(np.linalg.cond(pca_loader.reg_covmat_reps_mean))

            dof.append(pca_loader.n_comp)

        fig, ax1 = plotutils.subplots()
        ax1.plot(evr_range, l2_cond, label="L2 Condition Number")
        ax1.set_title(f"Dataset: {str(ds)}")
        ax1.set_xlabel("Explained Variance Ratio")
        ax1.set_ylabel("Covariance Matrix Condition Number")
        ax1.tick_params('y')

        # Draw horizontal line for threshold L2 condition number
        ax1.axhline(y=100, color='g', linestyle='--', label="Threshold L2 condition number")

        ax2 = ax1.twinx()

        # Plot the second dataset on the right y-axis
        ax2.plot(evr_range, dof, color="r", label="Degrees of freedom")
        ax2.set_ylabel('Degrees of freedom')
        ax2.tick_params('y')

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines = lines1 + lines2
        labels = labels1 + labels2
        ax1.legend(lines, labels, loc='upper left')

        yield fig


@figure
def xi_delta_histogram(normalized_delta_bias_data, title, lambda_value, label_hist=None):
    """
    Plot histogram of normalized delta regularized with PCA.

    Parameters
    ----------
    normalized_delta_bias_data: tuple

    label_hist: str
        summary description of multiclosure

    Returns
    -------
    fig
        Figure object
    """
    fig, ax = plotutils.subplots()
    size = np.shape(normalized_delta_bias_data[0])[0]
    ax.hist(normalized_delta_bias_data[0], density=True, bins=int(np.sqrt(size)), label=label_hist)
    ax.set_title(title + r", $\lambda:$" + f"{lambda_value}")
    ax.set_xlabel(r"$\delta$")
    ax.set_ylabel("Density")
    x = np.linspace(-3, 3, 100)
    y = norm.pdf(x, loc=0, scale=1)
    ax.plot(x, y, label="Standard gaussian")
    ax.legend()
    return fig


@table
def table_xi_indicator_function_data(bootstrapped_indicator_function_data):
    """
    Computes the bootstrap average and std of the indicator function for the data.

    Parameters
    ----------
    bootstrapped_indicator_function_data: tuple


    Returns
    -------
    pd.DataFrame
        DataFrame containing the average and std of the indicator function for the data.
    """
    indicator_list, mean_dof = bootstrapped_indicator_function_data

    # average over data and fits within the bootstrap samples
    mean_boot_vals = np.array([np.mean(boot_val) for boot_val in indicator_list])

    # take bootstrap expectation and variance
    mean_xi = np.mean(mean_boot_vals)
    std_xi = np.std(mean_boot_vals)

    records = [dict(data="full data", mean_dof=mean_dof, mean_xi=mean_xi, std_xi=std_xi)]

    df = pd.DataFrame.from_records(
        records, index="data", columns=("data", "mean_dof", "mean_xi", "std_xi")
    )
    return df
