"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""
import pandas as pd
import numpy as np
from scipy import stats

from reportengine.figure import figure, figuregen
from reportengine.table import table

from validphys import plotutils
from validphys.closuretest.inconsistent_closuretest.multiclosure_inconsistent import principal_components_dataset


@figuregen
def plot_bias_distribution_datasets(principal_components_bias_variance_datasets, each_dataset):
    """
    TODO
    """
    for pc_bias_var_dataset, ds in zip(principal_components_bias_variance_datasets, each_dataset):
        biases, variances, n_comp = pc_bias_var_dataset
        
        try:
            sqrt_rbv = np.sqrt(np.mean(biases) / np.mean(variances))
            vals = np.linspace(0, 3 * n_comp, 100)
            chi2_pdf = stats.chi2.pdf(vals, df=n_comp)
            chi2_cdf = lambda x: stats.chi2.cdf(x, df=n_comp)
            pvalue_ks = stats.kstest(biases, chi2_cdf).pvalue
            fig, ax = plotutils.subplots()
            ax.hist(biases, density=True, bins='auto', label=f"Bias {str(ds)}, p={pvalue_ks:.3f}")
            ax.plot(vals, chi2_pdf, label=f"chi2, dof={n_comp}")
            ax.plot([], [], 'ro', label=f"sqrt(Rbv) = {sqrt_rbv:.2f}")

            ax.legend()

            yield fig
        
        except:
            fig, ax = plotutils.subplots()
            ax.plot([], [], label=f"Dataset: {str(ds)}")
            ax.legend()
            yield fig

@figuregen
def plot_variance_distribution_datasets(
    principal_components_variance_distribution_datasets, each_dataset
):
    """
    TODO
    """

    for pc_var_dataset, ds in zip(
        principal_components_variance_distribution_datasets, each_dataset
    ):
        variances, n_comp = pc_var_dataset
        try:
            vals = np.linspace(0, 3 * n_comp, 100)
            chi2_pdf = stats.chi2.pdf(vals, df=n_comp)
            chi2_cdf = lambda x: stats.chi2.cdf(x, df=n_comp)

            for i, var_fit in enumerate(variances):
                pvalue_ks = stats.kstest(var_fit[0], chi2_cdf).pvalue

                fig, ax = plotutils.subplots()
                ax.hist(
                    var_fit[0],
                    density=True,
                    bins='auto',
                    label=f"Variance {str(ds)}, p={pvalue_ks:.3f}",
                )
                ax.plot(vals, chi2_pdf, label=f"chi2, dof={n_comp}")
                ax.set_title(f"Fit {i}")
                ax.legend()

                yield fig
        except:
            fig, ax = plotutils.subplots()
            yield fig

@table
def table_bias_variance_datasets(principal_components_bias_variance_datasets, each_dataset):
    """
    TODO
    """
    records = []
    for pc_bias_var_dataset, ds in zip(principal_components_bias_variance_datasets, each_dataset):
        biases, variances, n_comp = pc_bias_var_dataset
        
        try:
            bias = np.mean(biases)
            variance = np.mean(variances)
            rbv = bias / variance
            sqrt_rbv = np.sqrt(bias / variance)
            records.append(
            dict(
                dataset=str(ds),
                dof=n_comp,
                bias=bias,
                variance=variance,
                ratio=rbv,
                ratio_sqrt=sqrt_rbv,
            )
        )
        
        except:
            records.append(
            dict(
                dataset=str(ds),
                dof=n_comp,
                bias=bias,
                variance=variance,
                ratio=np.nan,
                ratio_sqrt=np.nan,
            ))
    
        
            

    df = pd.DataFrame.from_records(
            records,
            index="dataset",
            columns=("dataset", "dof", "bias", "variance", "ratio", "ratio_sqrt"),
        )
    df.columns = ["dof", "bias", "variance", "ratio", "sqrt(ratio)"]
    return df


@table
def datasets_bias_variance_pca_table(
    expected_datasets_fits_bias_variance_samples_pca, each_dataset
):
    """For each dataset calculate the expected bias and expected variance computed by
    first projecting the covariance matrix into a lower dimensional space trough PCA.

    """
    records = []
    for ds, (bias, var, ndata) in zip(
        each_dataset, expected_datasets_fits_bias_variance_samples_pca
    ):
        records.append(
            dict(
                dataset=str(ds),
                ndata=ndata,
                bias=bias,
                variance=var,
                ratio=bias / var,
                ratio_sqrt=np.sqrt(bias / var),
            )
        )
    df = pd.DataFrame.from_records(
        records,
        index="dataset",
        columns=("dataset", "ndata", "bias", "variance", "ratio", "ratio_sqrt"),
    )
    df.columns = ["ndata", "bias", "variance", "ratio", "sqrt(ratio)"]
    return df


@figure
def plot_sqrt_ratio_distribution_pca(dataset_fits_ratio_bias_variance_samples_pca):
    """"""
    sqrt_ratios = dataset_fits_ratio_bias_variance_samples_pca
    fig, ax = plotutils.subplots()
    ax.hist(sqrt_ratios, bins='auto', density=True, alpha=0.5, label="sqrt_ratio")
    ax.legend()
    return fig


@figure
def plot_variance_distribution_multi(multi_dataset_fits_bias_replicas_variance_samples, dataspecs):
    """
    histogram of the distribution of the variances (k) defined as the
    distance between central replica and replica (k) in units of the
    experimental covariance matrix.
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (_, variance_dist, _), spec in zip(
        multi_dataset_fits_bias_replicas_variance_samples, dataspecs
    ):
        label = spec['speclabel']

        ax.hist(variance_dist, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_variance_distribution_pca_multi(multi_dataset_fits_bias_variance_samples_pca, dataspecs):
    """
    like `plot_variance_distribution_multi`, but for variance computed with the covariance matrix
    predicted by the PDFs (and reduced with PCA).
    """
    return plot_variance_distribution_multi(multi_dataset_fits_bias_variance_samples_pca, dataspecs)


@figure
def plot_bias_distribution_multi(multi_dataset_fits_bias_replicas_variance_samples, dataspecs):
    """
    histogram of the distribution of the biases (l) defined as the
    distance between central replica and underlying law in units of the
    experimental covariance matrix.
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias_dist, _, _), spec in zip(
        multi_dataset_fits_bias_replicas_variance_samples, dataspecs
    ):
        label = spec['speclabel']

        ax.hist(bias_dist, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_bias_distribution_pca_multi(multi_dataset_fits_bias_variance_samples_pca, dataspecs):
    """
    like `plot_bias_distribution_multi`, but for variance computed with the covariance matrix
    predicted by the PDFs (and reduced with PCA).
    """
    return plot_bias_distribution_multi(multi_dataset_fits_bias_variance_samples_pca, dataspecs)


@figure
def plot_sqrt_ratio_distribution_pca(multi_dataset_fits_bias_variance_samples_pca, dataspecs):
    """
    histogram of the distribution of the sqrt ratio of bias and variance computed with
    the estimated "PDF" covariance matrix (reduced with PCA).
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    ratios = []
    labels = []
    for (bias_dist, variance_dist, _), spec in zip(
        multi_dataset_fits_bias_variance_samples_pca, dataspecs
    ):
        labels.append(spec['speclabel'])
        ratios.append(bias_dist / variance_dist)

    ax.hist(ratios, bins='auto', density=True, label=labels)
    ax.legend()
    return fig


@figure
def plot_multi_bias_distribution_bootstrap(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
    """
    histogram of the distribution of the biases obtained by bootstrapping with
    `fits_bootstrap_dataset_bias_variance` over the existing fits.

    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias, _), spec in zip(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
        label = spec['speclabel']

        ax.hist(bias, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_multi_ratio_bias_variance_distribution_bootstrap(
    multi_fits_bootstrap_dataset_bias_variance, dataspecs
):
    """
    histogram of the ratio bias variance distribution obtained by bootstrapping bias
    and variance with `fits_bootstrap_dataset_bias_variance`.

    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias, variance), spec in zip(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
        ratio = bias / variance
        label = spec['speclabel']

        ax.hist(ratio, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig

@figuregen
def plot_l2_condition_number(each_dataset, fits_pdf, variancepdf):
    """
    TODO
    """

    evr_range = np.linspace(0.90, 0.995, 10) 
    

    for ds in each_dataset:
        l2_cond = []
        dof = []

        for evr in evr_range:
            _, pc_reps, n_comp = principal_components_dataset(ds, fits_pdf, variancepdf, explained_variance_ratio=evr)
            
            covmat_pdf = np.cov(pc_reps)
            
            # compute condition number
            l2_cond.append(np.linalg.cond(covmat_pdf))
            dof.append(n_comp)
            


        fig, ax1 = plotutils.subplots(figsize=(15,4))
        ax1.plot(evr_range, l2_cond, "b-o", label="condition number")
        ax1.set_title(f"Dataset: {str(ds)}")
        ax1.set_xlabel("EVR")
        ax1.set_ylabel("Covariance Matrix Condition Number", color="b")
        ax1.tick_params('y', color="b")

        ax2 = ax1.twinx()

        # Plot the second dataset on the right y-axis
        ax2.plot(evr_range, dof, 'r-o', label="dof")
        ax2.set_ylabel('Degrees of freedom', color="r")
        ax2.tick_params('y', color="r")
        # ax1.legend()
        # ax2.legend()
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines = lines1 + lines2
        labels = labels1 + labels2
        ax1.legend(lines, labels, loc='upper left')


        yield fig