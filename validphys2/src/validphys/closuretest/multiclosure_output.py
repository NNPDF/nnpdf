"""
multiclosure_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for multiclosure estimators in the space of
data.

"""

import numpy as np
import pandas as pd
import scipy.special as special
from scipy.stats import norm

from reportengine.figure import figure, figuregen
from reportengine.table import table
from validphys import plotutils


@table
def table_bias_datasets(bias_datasets, each_dataset):
    """
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


@table
def bootstrapped_table_bias_datasets(bootstrapped_bias_datasets):
    """
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


@table
def bootstrapped_table_bias_data(bootstrapped_bias_data):
    """
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


@figuregen
def plot_xq2_data_prcs_maps(xq2_data_map, each_dataset):
    """
    Heat map of the ratio bias variance and xi quantile estimator for each datapoint
    in each dataset.

    Parameters
    ----------
    xq2_data_map: dictionary
    dictionary containing:
        - x coordinate
        - Q**2 coordinate
        - Ratio bias-variance
        - xi

    each_dataset: list

    Yields
    ------
    figure

    """
    keys = ["R_bv", "xi"]
    for j, elem in enumerate(xq2_data_map):

        for k in keys:
            if k == "R_bv":
                title = r"$R_{bv}$"
            if k == "xi":
                title = r"$\xi$"
            fig, ax = plotutils.subplots()
            im = ax.scatter(
                elem['x_coords'], elem['Q_coords'], c=(np.asarray(elem[k])), cmap='viridis'
            )
            fig.colorbar(im, label=title)
            ax.set_xscale('log')  # Set x-axis to log scale
            ax.set_yscale('log')  # Set y-axis to log scale
            ax.set_xlabel('x')
            ax.set_ylabel('Q2')
            ax.set_title(each_dataset[j].commondata.metadata.plotting.dataset_label)
            yield fig
