"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""

import pandas as pd
import numpy as np

from reportengine.figure import figuregen
from reportengine.table import table
from reportengine import collect

from validphys import plotutils
from validphys.closuretest.inconsistent_closuretest.multiclosure_inconsistent import (
    internal_multiclosure_dataset_loader_pca,
)


@table
def table_bias_variance_datasets(principal_components_bias_variance_datasets, principal_components_bias_variance_data, each_dataset):
    """
    Compute the bias, variance, ratio and sqrt(ratio) for each dataset
    and return a DataFrame with the results.

    Parameters
    ----------

    principal_components_bias_variance_datasets: list
        List of tuples containing the values of bias, variance and number of degrees of freedom
    
    principal_components_bias_variance_data: list
        Same of principal_components_bias_variance_datasets but for all the data

    each_dataset: list
        List of validphys.core.DataSetSpec

    Returns
    -------
    pd.DataFrame
        DataFrame containing the bias, variance, ratio and sqrt(ratio) for each dataset
    """
    records = []
    
    # First let's do the total
    biases_tot, variances_tot, n_comp_tot = principal_components_bias_variance_data
    bias_tot = np.mean(biases_tot)
    variance_tot = np.mean(variances_tot)
    rbv_tot = bias_tot / variance_tot
    # use gaussian uncertainty propagation
    delta_rbv_tot = np.sqrt(
            ((1 / variance_tot) * np.std(biases_tot)) ** 2 + (bias_tot / variance_tot**2 * np.std(variances_tot)) ** 2
        )
    sqrt_rbv_tot = np.sqrt(rbv_tot)
    delta_sqrt_rbv_tot = 0.5 * delta_rbv_tot / np.sqrt(rbv_tot)
    records.append(
            dict(
                dataset="Total",
                dof=n_comp_tot,
                bias=bias_tot,
                variance=variance_tot,
                ratio=rbv_tot,
                error_ratio=delta_rbv_tot,
                ratio_sqrt=sqrt_rbv_tot,
                error_ratio_sqrt=delta_sqrt_rbv_tot,
            )
        )
    # Now we do dataset per dataset
    for pc_bias_var_dataset, ds in zip(principal_components_bias_variance_datasets, each_dataset):
        biases, variances, n_comp = pc_bias_var_dataset
        bias = np.mean(biases)
        variance = np.mean(variances)
        rbv = bias / variance

        # use gaussian uncertainty propagation
        delta_rbv = np.sqrt(
            ((1 / variance) * np.std(biases)) ** 2 + (bias / variance**2 * np.std(variances)) ** 2
        )
        sqrt_rbv = np.sqrt(bias / variance)
        delta_sqrt_rbv = 0.5 * delta_rbv / np.sqrt(rbv)

        records.append(
            dict(
                dataset=str(ds),
                dof=n_comp,
                bias=bias,
                variance=variance,
                ratio=rbv,
                error_ratio=delta_rbv,
                ratio_sqrt=sqrt_rbv,
                error_ratio_sqrt=delta_sqrt_rbv,
            )
        )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=(
            "dataset",
            "dof",
            "bias",
            "variance",
            "ratio",
            "error_ratio",
            "ratio_sqrt",
            "error_ratio_sqrt",
        ),
    )
    df.columns = [
        "dof",
        "bias",
        "variance",
        "ratio",
        "error ratio",
        "sqrt(ratio)",
        "error sqrt(ratio)",
    ]

    return df


lambdavalues_table_bias_variance_datasets = collect(
    "table_bias_variance_datasets", ("lambdavalues",)
)


@table
def bootstrapped_table_bias_variance_datasets(bootstrapped_principal_components_bias_variance_data):
    """
    Compute the bias, variance, ratio and sqrt(ratio) for each dataset
    and return a DataFrame with the results.
    Uncertainty on ratio and sqrt ratio is computed by Gaussian error propagation
    of the bootstrap uncertainty on bias and variance.
    """
    records = []
    for boot_ds in bootstrapped_principal_components_bias_variance_data:
        df = boot_ds
        mean_bias = df["bias"].mean()
        mean_variance = df["variance"].mean()
        mean_ratio = mean_bias / mean_variance
        bootstrap_unc_bias = df["bias"].std()

        # gaussian error propagation for the ratio of the means uncertainty
        # only consider bias as source of uncertainty for the ratio (variance is almost constant)
        bootstrap_unc_ratio = np.sqrt((1 / mean_variance * bootstrap_unc_bias) ** 2)
        sqrt_ratio = np.sqrt(mean_ratio)
        # gaussian error propagation for the sqrt of the ratio
        bootstrap_unc_sqrt_ratio = 0.5 * bootstrap_unc_ratio / np.sqrt(mean_ratio)

        records.append(
            dict(
                dataset=df["dataset"].iloc[0],
                mean_dof=df.n_comp.mean(),
                bias=mean_bias,
                variance=mean_variance,
                ratio=mean_ratio,
                error_ratio=bootstrap_unc_ratio,
                ratio_sqrt=sqrt_ratio,
                error_ratio_sqrt=bootstrap_unc_sqrt_ratio,
            )
        )

    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=(
            "dataset",
            "mean_dof",
            "bias",
            "variance",
            "ratio",
            "error_ratio",
            "ratio_sqrt",
            "error_ratio_sqrt",
        ),
    )
    df.columns = [
        "mean_dof",
        "bias",
        "variance",
        "ratio",
        "error_ratio",
        "ratio_sqrt",
        "error_ratio_sqrt",
    ]

    return df


@figuregen
def plot_lambdavalues_bias_variance_values(
    lambdavalues_table_bias_variance_datasets, lambdavalues, each_dataset
):
    """
    Plot sqrt of ratio bias variance as a function of lambda for each dataset.

    Parameters
    ----------
    lambdavalues_table_bias_variance_datasets: list
        list of data frames computed as per table_bias_variance_datasets.

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
            df = lambdavalues_table_bias_variance_datasets[i]
            df = df[df.index == str(ds)]

            ax.errorbar(
                lambdavalue["lambda_value"],
                df["sqrt(ratio)"].values,
                yerr=df["error sqrt(ratio)"].values,
                color="blue",
                fmt='o',
            )
            ax.hlines(1, xmin=0, xmax=1.0, color="red", linestyle="--")
            ax.set_ylabel(r"$R_{bv}$")
            ax.set_xlabel(r"$\lambda$")

        ax.set_title(f"{ds.commondata.metadata.plotting.dataset_label}")

        yield fig


internal_multiclosure_data_collected_loader = collect(
    "internal_multiclosure_dataset_loader", ("data",)
)


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

    internal_multiclosure_data_loader: list
        list of internal_multiclosure_dataset_loader objects


    Yields
    ------
    fig
        Figure object
    """

    # Explained variance ratio range
    evr_range = np.linspace(evr_min, evr_max, evr_n)

    for internal_multiclosure_dataset_loader, ds in zip(
        internal_multiclosure_data_collected_loader, each_dataset
    ):
        l2_cond = []
        dof = []

        for evr in evr_range:

            pca_loader = internal_multiclosure_dataset_loader_pca(
                internal_multiclosure_dataset_loader, explained_variance_ratio=evr
            )

            # if the number of principal components is 1 then the covariance matrix is a scalar
            # and the condition number is not computed
            if pca_loader.n_comp == 1:
                l2_cond.append(np.nan)
            else:
                l2_cond.append(np.linalg.cond(pca_loader.covmat_pca))

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
