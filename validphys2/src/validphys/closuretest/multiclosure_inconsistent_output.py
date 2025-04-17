"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""

import numpy as np

from reportengine.figure import figuregen, figure
from reportengine import collect

from validphys import plotutils
from validphys.closuretest.multiclosure import regularized_multiclosure_dataset_loader


"""
Collects bootstrapped_table_bias_datasets over multiple lambda values dataspecs.
"""
lambdavalues_bootstrapped_table_bias_datasets = collect(
    "bootstrapped_table_bias_datasets", ("lambdavalues",)
)


"""
Collects bootstrapped_table_bias_data over multiple lambda values dataspecs.
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
