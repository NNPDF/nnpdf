"""
plots.py

Plots for the paramfits package.
"""
import logging
import numbers

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from reportengine.floatformatting import format_number
from reportengine.checks import make_argcheck, CheckError, check_positive
from reportengine.figure import figure, figuregen

from validphys.plotutils import plot_horizontal_errorbars, barplot, kde_plot, marker_iter_plot
from validphys.paramfits.dataops import check_dataset_items, get_parabola

log = logging.getLogger(__name__)

@figure
def plot_fits_as_profile(fits_as, fits_total_chi2, suptitle=None):
    """Plot the total central chi² as a function of the value of α_s.
    Note that this plots as a function of the key "AlphaS_MZ" in the LHAPDF
    file, which is annoyingly *not* α_s(MZ) for Nf<5."""
    fig, ax = plt.subplots()
    alphas = fits_as
    #Could be a transposed data frame
    fits_total_chi2 = np.ravel(fits_total_chi2)
    ax.plot(alphas, fits_total_chi2)
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²$')
    if suptitle:
        fig.suptitle(suptitle)
    return fig

@figure
def plot_as_central_parabola(
        fits_as,
        as_central_parabola,
        suptitle, ndata,
        parabolic_as_determination_for_total,
        markmin:bool=False):
    """Plot a parabola with the central chi² per number of points, marking
    the chi² at the total best fit."""
    fig,ax = plt.subplots()
    asarr = np.linspace(min(fits_as), max(fits_as), 100)
    as_central_parabola = np.asarray(as_central_parabola)
    ndata = int(ndata)
    ax.plot(asarr, np.polyval(as_central_parabola, asarr)/ndata)

    best_as = parabolic_as_determination_for_total[0].location
    chi2_at_best = np.polyval(as_central_parabola, best_as)/ndata
    ax.scatter(best_as, chi2_at_best)
    ax.annotate(format_number(chi2_at_best, 3), (best_as, chi2_at_best))
    if markmin:
        ax.axvline(-as_central_parabola[1]/2/as_central_parabola[0])
    ax.set_ylabel(r'$\chi^2/N_{data}$')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_title(f"{suptitle}")
    return fig

@figure
def plot_as_cummulative_central_chi2(fits_as,
                                     as_datasets_central_parabolas,
                                     central_by_dataset_suptitle):
    """Plot the cumulative total chi² for each of the datasets"""
    fig,ax = plt.subplots()
    nx = 100
    asarr = np.linspace(min(fits_as), max(fits_as), nx)
    last = np.zeros(nx)
    for (p, label) in zip(as_datasets_central_parabolas, central_by_dataset_suptitle):
        val = last + np.polyval(p, asarr)
        ax.fill_between(asarr, last, val, label=label)
        last = val
    ax.legend()
    ax.set_ylabel(r'$\chi^2$')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylim(0)
    ax.set_xlim(asarr[[0,-1]])

    return fig

@figure
def plot_as_cummulative_central_chi2_diff(fits_as,
                                     as_datasets_central_parabolas,
                                     central_by_dataset_suptitle,
                                     parabolic_as_determination_for_total,
                                     ymax:(numbers.Real, type(None))=None):
    """Plot the cumulative difference between the χ² at the best global
    αs fit and the χ² at αs. If the difference is negative, it is set to zero.
    """
    fig,ax = plt.subplots()
    nx = 100
    best_as = parabolic_as_determination_for_total[0].location
    asarr = np.linspace(min(fits_as), max(fits_as), nx)
    last = np.zeros(nx)
    for (p, label) in zip(as_datasets_central_parabolas, central_by_dataset_suptitle):
        delta = np.polyval(p, asarr) - np.polyval(p, best_as)
        delta[delta<0] = 0
        val = last + delta
        ax.fill_between(asarr, last, val, label=label)
        last = val
    ax.legend()
    ax.set_ylabel(r'$\chi^2 - \chi^2_{{\rm best}\ \alpha_S}$')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylim(0)
    if ymax is not None:
        ax.set_ylim(ymax=ymax)
    ax.set_xlim(asarr[[0,-1]])

    return fig

@figure
def plot_as_cummulative_central_chi2_diff_underflow(
                                     fits_as,
                                     as_datasets_central_parabolas,
                                     central_by_dataset_suptitle,
                                     parabolic_as_determination_for_total,
                                     ymax:(numbers.Real, type(None))=None):
    """Plot the cumulative difference between the χ² at the best global
    αs fit and the χ² at αs. If the difference is negative, it is set to zero.
    """
    fig,ax = plt.subplots()
    nx = 100
    best_as = parabolic_as_determination_for_total[0].location
    asarr = np.linspace(min(fits_as), max(fits_as), nx)

    #ordering = [np.polyval(p, best_as + 0.001) - np.polyval(p, best_as) for p in as_datasets_central_parabolas]
    #ordered_indexes = sorted(range(len(ordering)), key=ordering.__getitem__)
    #as_datasets_central_parabolas = [as_datasets_central_parabolas[i] for i in ordered_indexes]
    #central_by_dataset_suptitle = [central_by_dataset_suptitle[i] for i in ordered_indexes]



    last = np.zeros(nx)
    for (p, label) in zip(as_datasets_central_parabolas, central_by_dataset_suptitle):
        delta = np.polyval(p, asarr) - np.polyval(p, best_as)
        delta[delta<0] = 0
        val = last + delta
        ax.fill_between(asarr, last, val, label=label)
        last = val
    last = np.zeros(nx)
    for (p, label) in zip(as_datasets_central_parabolas, central_by_dataset_suptitle):
        delta = np.polyval(p, asarr) - np.polyval(p, best_as)
        delta[delta>0] = 0
        val = last + delta
        ax.fill_between(asarr, last, val)
        last = val
    ax.legend()
    ax.set_ylabel(r'$\chi^2 - \chi^2_{{\rm best}\ \alpha_S}$')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_xlim(asarr[[0,-1]])
    if ymax is not None:
        ax.set_ylim(ymax=ymax)
    ymin = np.min(last)
    ax.set_ylim(ymin=ymin)

    return fig

@figure
def plot_as_cummulative_central_chi2_diff_negative(
                                     fits_as,
                                     as_datasets_central_parabolas,
                                     central_by_dataset_suptitle,
                                     parabolic_as_determination_for_total):
    """Plot the cumulative difference between the χ² at the best global
    αs fit and the χ² at αs. If the difference is negative, it is set to zero.
    """
    """Plot the cumulative difference between the χ² at the best global
    αs fit and the χ² at αs. If the difference is negative, it is set to zero.
    """
    fig,ax = plt.subplots()
    nx = 100
    best_as = parabolic_as_determination_for_total[0].location
    asarr = np.linspace(min(fits_as), max(fits_as), nx)
    last = np.zeros(nx)
    for (p, label) in zip(as_datasets_central_parabolas, central_by_dataset_suptitle):
        delta = np.polyval(p, asarr) - np.polyval(p, best_as)
        delta[delta>0] = 0
        val = last + delta
        ax.fill_between(asarr, last, val, label=label)
        last = val
    #Seems that fill_between doesn't compute for the legend location
    ax.legend(loc='lower left')
    ax.set_ylabel(r'$\chi^2 - \chi^2_{{\rm best}\ \alpha_S}$')
    ax.set_xlabel(r'$\alpha_S$')
    ymin = np.min(last)
    ax.set_ylim(ymin=ymin)
    ax.set_xlim(asarr[[0,-1]])

    return fig
@figuregen
def plot_dataspecs_central_parabolas(
        dataspecs_as_central_parabolas_map,
        dataspecs_fits_as):
    """Plot the parabolas resulting from the chi² of the mean PDF to the data,
    as a function of alpha_S. Yield one plot per ``dataset_item`` comparing
    several dataspecs."""
    #Note that cannot cast as a matrix if shapes are different
    limits = min(map(min, dataspecs_fits_as)), max(map(max, dataspecs_fits_as))
    alphasaxis = np.linspace(*limits, 100)
    for dsname, dataspecs_parabolas in dataspecs_as_central_parabolas_map.items():
        fig, ax = plt.subplots()
        for label, parabola in dataspecs_parabolas.items():
            ax.plot(alphasaxis, np.polyval(parabola, alphasaxis), label=label)

        ax.set_ylabel(r'$\chi^2/N_{dat}$')
        ax.set_xlabel(r'$\alpha_S$')
        #ax.set_ylim(0)
        ax.set_title(dsname)
        ax.set_xlim(alphasaxis[[0,-1]])

        plt.legend()

        yield fig

@figure
def plot_as_datasets_pseudoreplicas_chi2(as_datasets_pseudoreplicas_chi2):
    """Plot the error bars of the αs determination from pseudoreplicas
    by dataset item. Note that this only has meaning of preferred
    value for "Total", and the rest of the values are the minima of
    the partial χ²."""
    data, names = zip(*as_datasets_pseudoreplicas_chi2)
    cv, err = zip(*[(np.mean(dt), np.std(dt)) for dt in data])
    fig, ax = plot_horizontal_errorbars([cv], [err], names)
    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ from pseudoreplicas")
    return fig

@figure
def plot_as_exepriments_central_chi2(as_datasets_central_chi2):
    """Plot the error bars of the αs determination from central χ²
    by dataset item. Note that this only has meaning of preferred
    value for "Total", and the rest of the values are the minima of
    the partial χ²."""
    data, names = zip(*as_datasets_central_chi2)
    cv, err = zip(*data)
    fig, ax = plot_horizontal_errorbars([cv], [err], names)
    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ from central chi²")
    return fig


@figure
def plot_as_datasets_compare(as_datasets_pseudoreplicas_chi2,
                             as_datasets_central_chi2,
                             marktotal:bool=True):
    """Plot the result of ``plot_as_datasets_pseudoreplicas_chi2`` and
    ``plot_as_exepriments_central_chi2`` together."""
    datapseudo, namespseudo = zip(*as_datasets_pseudoreplicas_chi2)
    cvpseudo, errpseudo = zip(*[(np.mean(dt), np.std(dt)) for dt in datapseudo])


    datacentral, namescentral = zip(*as_datasets_central_chi2)
    cvcentral, errcentral = zip(*datacentral)

    if namespseudo != namescentral:
        raise RuntimeError("Names do not coincide")

    fig, ax = plot_horizontal_errorbars(
        [cvcentral, cvpseudo], [errcentral, errpseudo], namescentral,
        [r'Central $\chi^2$', r'Pseudoreplica $\chi^2$']
    )
    if marktotal:
        try:
            pos = namespseudo.index('Total')
        except ValueError:
            log.error("Asked to mark total, but it was not provided.")
        else:
            ax.axvline(cvcentral[pos], color='C0', linewidth=0.5, linestyle='--')
            ax.axvline(cvpseudo[pos], color='C1', linewidth=0.5, linestyle='--')

    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ determination")
    ax.legend()
    return fig

@figure
def plot_dataspecs_as_value_error(dataspecs_as_value_error_table_impl,
        dataspecs_fits_as,
        marktotal:bool=True, fix_limits:bool=True):
    """
    Plot the result for each dataspec of the pseudoreplica alpha_s
    determination based on the partial chi² for each ``dataset_item``.

    If ``marktotal`` is True, a vertical line will appear marking the position
    of the best fit.

    If ``fix_limits`` is True, the limits of the plot will span all the fitted
    values. Otherwise an heuristic will be used.

    """

    df = dataspecs_as_value_error_table_impl
    datalabels = df.columns.levels[0]
    catlabels = list(df.index)
    cvs = df.loc[:, (slice(None), 'mean')].T.values
    errors = df.loc[:, (slice(None), 'error')].T.values

    if fix_limits:
        minlim = min(min(x for x in dataspecs_fits_as))
        maxlim = max(max(x for x in dataspecs_fits_as))
        lims = minlim, maxlim
    else:
        lims = None


    fig, ax = plot_horizontal_errorbars(
        cvs, errors, catlabels,
        datalabels,
        xlim = lims

    )

    if marktotal:
        try:
            pos = catlabels.index('Total')
        except ValueError:
            log.error("Asked to mark total, but it was not provided.")
        else:
            for i,cv in enumerate(cvs):
                ax.axvline(cv[pos], color=f'C{i%10}', linewidth=0.5, linestyle='--')

    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ determination")
    ax.legend()
    return fig

@figure
def plot_dataspecs_as_value_error_comparing_with_central(
        dataspecs_as_value_error_table_impl,
        as_datasets_central_chi2,
        dataspecs_fits_as,
        speclabel,
        marktotal:bool=True, fix_limits:bool=True):
    """
    This is an aberration we need to do for the paper plots. It compares the
    central (old) and new partial chi².
    """

    df = dataspecs_as_value_error_table_impl
    catlabels = list(df.index)
    replica_cvs = df.loc[:, (speclabel, 'mean')].T.values
    replica_errors = df.loc[:, (speclabel, 'error')].T.values

    datacentral, namescentral = zip(*as_datasets_central_chi2)
    cvcentral, errcentral = zip(*datacentral)

    if fix_limits:
        minlim = min(min(x for x in dataspecs_fits_as))
        maxlim = max(max(x for x in dataspecs_fits_as))
        lims = minlim, maxlim
    else:
        lims = None


    cvs = np.c_[cvcentral, replica_cvs].T
    errors = np.c_[errcentral, replica_errors].T



    fig, ax = plot_horizontal_errorbars(
        cvs, errors, catlabels,
        [r'$\Delta \chi^2=1$', 'Replicas'],
        xlim = lims

    )

    if marktotal:
        try:
            pos = catlabels.index('Total')
        except ValueError:
            log.error("Asked to mark total, but it was not provided.")
        else:
            for i,cv in enumerate(cvs):
                ax.axvline(cv[pos], color=f'C{i%10}', linewidth=0.5, linestyle='--')

    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ determination")
    ax.legend()
    return fig

@make_argcheck
def _check_first_is_total(fits_central_chi2_by_experiment_and_dataset):
    l = fits_central_chi2_by_experiment_and_dataset
    if not l or l[0]['experiment_label'] != 'Total':
        raise CheckError("Expecting that the first experiment is the total. You may need to set prepend_total=True")

@figure
@_check_first_is_total
def plot_as_value_error_central(as_datasets_central_chi2,
         marktotal:bool=True):
    """Plot the result of ``plot_as_datasets_pseudoreplicas_chi2`` and
    ``plot_as_exepriments_central_chi2`` together."""

    datacentral, namescentral = zip(*as_datasets_central_chi2)
    cvcentral, errcentral = zip(*datacentral)

    fig, ax = plot_horizontal_errorbars(
         [cvcentral], [errcentral], namescentral,
         [r'Central']
    )
    ax.axvline(cvcentral[0], color=f'C{0}', linewidth=0.5, linestyle='--')

    ax.set_xlabel(r"$\alpha_S$")
    ax.set_xlim(0.110,0.125)
    ax.set_title(r"$\alpha_S$ determination")
    ax.legend()
    return fig

# Pull plots
def _pulls_func(cv,alphas_global,error,error_global):
    """Small definition to compute pulls"""
    return error_global*((cv-alphas_global)/(error**2))


@figure
@_check_first_is_total
def plot_pulls_central(as_datasets_central_chi2,hide_total:bool=True):
    """ Plots the pulls per experiment for the central results """

    data, names = zip(*as_datasets_central_chi2)
    cv, err = zip(*data)
    pulls = list()
    if hide_total:
        for i in range(1,len(cv)):
            pulls.append(_pulls_func(cv[i],cv[0],err[i],err[0]))
            names = [x for x in names if x!='Total']
    else:
        for i in range(0,len(cv)):
            pulls.append(_pulls_func(cv[i],cv[0],err[i],err[0]))

    fig, ax = barplot(pulls, names, " ", orientation="horizontal")
    ax.legend()

    return fig

@figure
@_check_first_is_total
def plot_pull_gaussian_fit_central(as_datasets_central_chi2,
        dataspecs_fits_as,dataspecs_speclabel,hide_total:bool=True):

    """Bins the pulls and overlays
    the normalised gaussian fit and KDE to the histogram of pulls"""

    data, names = zip(*as_datasets_central_chi2)
    cv, err = zip(*data)
    pulls = list()

    if hide_total:
        for i in range(1,len(cv)):
            pulls.append(_pulls_func(cv[i],cv[0],err[i],err[0]))
            names = [x for x in names if x!='Total']
    else:
        for i in range(0,len(cv)):
            pulls.append(_pulls_func(cv[i],cv[0],err[i],err[0]))

    mean_pulls = np.mean(pulls)
    std_dev = np.std(pulls)
    x = np.linspace(min(pulls),max(pulls), 100)
    kde_pulls = stats.gaussian_kde(pulls, bw_method='silverman')
    fig, ax = plt.subplots()

    #ax.set_title(f"Histogram of pulls for {label} dataset")
    ax.set_xlabel(r"Pull")
    ax.plot(x, kde_pulls(x), label="Kernel Density Estimation of pulls")
    ax.hist(pulls,normed=True,bins=4)
    ax.grid(False)
    normdist = stats.norm(mean_pulls, std_dev)
    ax.plot(x, normdist.pdf(x),label="Normalised gaussian fit")
    ax.legend()

    return fig

@figure
def plot_pull_plots_global_min(dataspecs_as_value_error_table_impl,
        dataspecs_fits_as,dataspecs_speclabel,hide_total:bool=True):

    """Plots the pulls of individual experiments as a barplot."""

    df = dataspecs_as_value_error_table_impl
    tots_error = df.loc['Total', (slice(None), 'error')].values
    tots_mean = df.loc['Total', (slice(None), 'mean')].values

    if hide_total:
        df = df.loc[df.index != 'Total']

    catlabels = list(df.index)
    cvs = df.loc[:, (slice(None), 'mean')].values
    errors = df.loc[:, (slice(None), 'error')].values

    pulls = _pulls_func(cvs,tots_mean,errors,tots_error).T

    pulls = np.asarray(pulls,dtype=float)
    #Compute sums
    #sp = np.atleast_2d(np.nansum(pulls, axis=1)).T
    #pulls =np.concatenate([pulls, sp], axis=1)

    #fig, ax = barplot(pulls, [*catlabels, 'sum'], dataspecs_speclabel, orientation="horizontal")
    fig, ax = barplot(pulls, catlabels, dataspecs_speclabel, orientation="horizontal")


    #ax.set_title(r"Pulls per experiment")
    ax.legend()
    return fig

@make_argcheck
def _check_two_speclabels(dataspecs_speclabel):
    if len(dataspecs_speclabel) != 2:
        raise CheckError("Need 2 data specs")


@figure
@_check_two_speclabels
def alphas_shift(
    dataspecs_as_value_error_table_impl,
    dataspecs_quad_table_impl,
    dataspecs_ndata_table,
    dataspecs_dataset_ndata,
    dataspecs_fits_as,
    dataspecs_speclabel,
    hide_total:bool=True,
    ndata_weight:bool=False):

    """Plots NNLO - NLO alphas values for each experiment - i.e.
        the shift in the best fit alphas for each process (as it currently
        stands...) wrt the global best fit alphas at NLO or NNLO.
        Also contains some computations for estimating MHOU, using either
        the number of data points per experiment/process (ndata)
        or the quadratic coefficient of the parabolic fit (quad_weights)"""

    df1 = dataspecs_ndata_table
    df = dataspecs_as_value_error_table_impl
    df2 = dataspecs_quad_table_impl


    tots_mean = df.loc['Total', (slice(None), 'mean')].values

    if hide_total:
        df = df.loc[df.index != 'Total']
        df1 = df1.loc[df1.index != 'Total']
        df2 = df2.loc[df2.index != 'Total']


    cvs = df.loc[:, (slice(None), 'mean')].T.values
    quad_weights = df2.loc[:, (slice(None), 'mean')].T.values

    catlabels = list(df.index)

    alphas_shift = []
    nnlo_alphas_global_shift = []
    nlo_alphas_global_shift = []
    nnlo2_alphas_global_shift = []

    for i in range(0,len(cvs[0])):
        alphas_shift.append(cvs[1][i]-cvs[0][i])
        nnlo2_alphas_global_shift.append((cvs[1][i]-tots_mean[1])**2)
        nlo_alphas_global_shift.append((cvs[0][i]-tots_mean[0])**2)
        nnlo_alphas_global_shift.append(cvs[1][i]-tots_mean[1])

    weights_nlo = []
    weights_nnlo = []
    weights_nnlo_sq = []
    weights_nlo_sq = []

    if ndata_weight:

        ndataptsnlo = df1.iloc[:,0]
        ndataptsnnlo = df1.iloc[:,1]

        for i in range(0,len(ndataptsnlo)):
            weights_nlo.append(ndataptsnlo[i])
            weights_nnlo.append(ndataptsnnlo[i])

    else:
        for i in range(0,len(quad_weights.T)):
            weights_nlo.append(quad_weights[0][i])
            weights_nnlo.append(quad_weights[1][i])
            weights_nnlo_sq.append((quad_weights[1][i])**2)
            weights_nlo_sq.append((quad_weights[0][i])**2)


    term1, term2 = dataspecs_speclabel


    fig, ax = barplot(alphas_shift, catlabels, " ", orientation = "horizontal")
    ax.set_title(f"{term2} - {term1} shifts")
    ax.legend()
    return fig

@figuregen
def plot_pull_gaussian_fit_pseudo(dataspecs_as_value_error_table_impl,
        dataspecs_fits_as,dataspecs_speclabel,hide_total:bool=True):

    """Bins the pulls computed in pull_plots_global_min and overlays
    the normalised gaussian fit and KDE to the histogram of pulls"""

    df = dataspecs_as_value_error_table_impl
    tots_error = df.loc['Total', (slice(None), 'error')].T.values
    tots_mean = df.loc['Total', (slice(None), 'mean')].T.values

    if hide_total:
        df = df.loc[df.index != 'Total']

    cvs = df.loc[:, (slice(None), 'mean')].T.values
    errors = df.loc[:, (slice(None), 'error')].T.values

    for label, i in zip(dataspecs_speclabel, range(len(cvs))):
        pulls = _pulls_func(cvs[i],tots_mean[i],errors[i],tots_error[i])

        mean_pulls = np.mean(pulls)
        std_dev = np.std(pulls)
        x = np.linspace(min(pulls),max(pulls), 100)

        kde_pulls = stats.gaussian_kde(pulls, bw_method='silverman')
        fig, ax = plt.subplots()

        #ax.set_title(f"Histogram of pulls for {label} dataset")
        ax.set_xlabel(r"Pull")
        ax.plot(x, kde_pulls(x), label="Kernel Density Estimation of pulls")
        ax.hist(pulls,normed=True,bins=4)
        ax.grid(False)
        normdist = stats.norm(mean_pulls, std_dev)
        ax.plot(x, normdist.pdf(x) ,label="Normalised gaussian fit")
        ax.legend()

        yield fig

@figure
def plot_fitted_replicas_as_profiles_matched(fits_as,
        fits_replica_data_with_discarded_replicas,
        parabolic_as_determination, suptitle=None):
    """Plot chi²(as) keeping the replica nnfit index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """
    alphas = fits_as


    minimums = parabolic_as_determination

    table = fits_replica_data_with_discarded_replicas.values

    fig, ax = plt.subplots()

    from matplotlib.collections import LineCollection

    lc = LineCollection([list(x for x in zip(alphas, t) if np.isfinite(x[1])) for t in table])
    lc.set_array(minimums.data)
    lc.set_clim(*np.percentile(minimums.data, (5,95)))
    ax.add_collection(lc)
    ax.set_xlim(min(alphas), max(alphas))
    ax.set_ylim(np.nanmin(table), np.nanmax(table))
    fig.colorbar(lc, label=r"Preferred $\alpha_S$")
    ax.set_title(rf"$\alpha_S$ = ${minimums.location:.4f} \pm {minimums.scale:.4f}$ N={len(minimums.data)}")
    if suptitle:
        fig.suptitle(suptitle)

    #ax.plot(alphas, np.array(table).T, color='#ddddcc')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi^2$')
    return fig

@figuregen
@check_dataset_items
def plot_dataspecs_pseudoreplica_means(
        dataspecs_chi2_by_dataset_dict,
        dataspecs_speclabel,
        dataset_items:(list, type(None))=None,
        ):
    """ Plot the mean chi² from data to pseudoreplica, over replicas in a fit
    and comparing dataspecs.
    """
    if dataset_items is None:
        dataset_items = list(dataspecs_chi2_by_dataset_dict)


    for it in dataset_items:
        fig, ax = plt.subplots()
        for label, tb in zip(dataspecs_speclabel, dataspecs_chi2_by_dataset_dict[it]):
            if tb is None:
                #Advance color cycle so all dataspecs get the same color
                ax._get_lines.get_next_color()
                continue
            m = tb.mean(axis=0)
            ax.plot(np.asarray(m.index),np.asarray(m), label=label)
        ax.set_title(it)
        ax.set_xlabel(r'$\alpha_S$')
        ax.set_ylabel(r'$\left<\chi^2\right>$')
        ax.legend()

        yield fig



#TODO: This should take parabolic_as_determination with the as_transforms
#and so on, rather that refitting here.
@figuregen
@check_dataset_items
@check_positive('examples_per_item')
def plot_dataspecs_parabola_examples(
        dataspecs_chi2_by_dataset_dict,
        dataspecs_speclabel,
        dataset_items:(list, type(None))=None,
        examples_per_item:int = 2,
        random_seed:int = 0,
        ):
    """Sample ``examples_per_item`` replica_indexes for each of the
    ``dataset_items``. Yield a plot with the parabolic fit, as resolved for
    each of the dataspecs. The random state is local to the function and
    controlled by ``random_seed``."""
    random_state = np.random.RandomState(random_seed)
    if dataset_items is None:
        dataset_items = list(dataspecs_chi2_by_dataset_dict)


    for it in dataset_items:
        table = pd.concat(dataspecs_chi2_by_dataset_dict[it],
                         join='inner', axis=1, keys=dataspecs_speclabel)
        if examples_per_item > len(table):
            log.warning("Length of table is less than examples_per_item")
            sampled = table
        else:
            sampled = table.sample(examples_per_item, random_state=random_state)


        for index, row in sampled.iterrows():

            fig, ax = plt.subplots()
            im = marker_iter_plot()
            ax.set_title(f"Parabola example for {it} (nnfit_index={index})")
            for i, (label, vals) in enumerate(row.groupby(level=0)):
                asvals = vals.index.get_level_values(1)
                color = f'C{i}'
                y = vals.values
                ax.plot(asvals, y, **next(im), label=label,
                         color=color, linestyle='none', lw=0.5)
                a,b,c = parabola = get_parabola(asvals, y)
                ax.plot(asvals, np.polyval(parabola,asvals), color=color,
                         linestyle='--', lw=0.5)
                m = -b/2/a
                if asvals[0] < m < asvals[-1]:
                    ax.axvline(m ,  color=color, lw=0.4)
                ax.set_xlabel(r'$\alpha_S$')
                ax.set_ylabel(r'$\chi^2$')
            ax.legend(  )

            yield fig

@figure
def plot_as_distribution(parabolic_as_determination, suptitle):
    """Histograms of the values of alphas produced, with the datapoints in
    an array as sticks on an axis"""

    distribution = parabolic_as_determination.data

    fig, ax = plt.subplots()

    kde_plot(distribution)
    ax.legend()
    ax.set_title(f"{suptitle}")
    ax.set_xlabel(r"$\alpha_S$")
    return fig

@figure
def plot_total_as_distribution_dataspecs(
        dataspecs_parabolic_as_determination_for_total,
        dataspecs_speclabel,
        ):
    """Compare the total alpha_s distributions across dataspecs.
    See ``plot_as_distribution``."""
    fig, ax = plt.subplots()
    for dist, label in zip(
            dataspecs_parabolic_as_determination_for_total,
            dataspecs_speclabel):
        #Remember that *_for_total is a len 1 list, so take the first element.
        kde_plot(dist[0].data, ax=ax, label=label)
    ax.set_xlabel(r"$\alpha_S$")
    ax.legend()
    return fig

@figure
def plot_poly_as_fit(fits_as,
        fits_replica_data_correlated, max_ndiscarded:int=4, polorder:int=2,
        suptitle=None):
    """Plot a polynomial fit of chi²(as) of `degree polorder`, keeping the
    replica index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """

    alphas = fits_as
    df = fits_replica_data_correlated
    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded
    table = table[filt]
    table = table.values
    fig, ax = plt.subplots()

    minimums = []
    asarr = np.asarray(alphas)
    for i, row in enumerate(table):
        filt =  np.isfinite(row)
        fit = np.polyfit(asarr[filt], row[filt], polorder)
        extrema = np.roots(np.polyder(fit))

        #Cast away the zero complex part
        candidates = np.real(extrema[np.isreal(extrema)])
        #candidates = candidates[(candidates > alphas[0]) & (candidates<alphas[-1])]
        if len(candidates):
            minimums.append(candidates[np.argmin(np.polyval(fit, candidates))])
            color = f'C{i%10}'
            ax.plot(alphas, np.polyval(fit, alphas), color=color)
            ax.plot(asarr[filt], row[filt], 'o', color=color)
    minimums = np.asarray(minimums)
    ax.set_xlim(min(alphas), max(alphas))
    ax.set_ylim(np.nanmin(table), np.nanmax(table))
    plt.title(rf"$\alpha_S$ from order {polorder} polynomial fit "
              rf"= ${np.mean(minimums):.4f} \pm {np.std(minimums):.4f}$"
              rf"N={len(minimums)}")
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²/N_{dat}$')
    if suptitle:
        fig.suptitle(suptitle)
    return fig


#TODO: This should probably be done more granularly and not here

def _reduce_mean_parabola(df):
    if df is None:
        return None
    a,b,c = get_parabola(df.columns, df.mean(axis=0))
    if a < 0:
        raise ValueError("Expecting convex parabola")
    minimum = -b/2/a
    error  = 1/np.sqrt(a)
    return minimum, error

@figure
def plot_mean_pulls(dataspecs_chi2_by_dataset_dict, dataspecs_speclabel):
    """Compute the pulls from the sum of the parabolas."""
    d = dataspecs_chi2_by_dataset_dict
    dtotal = d['Total']
    er_val_global = [_reduce_mean_parabola(v) for v in dtotal]
    ks = set(d.keys())
    ks.remove('Total')
    l = [0]*len(dataspecs_speclabel)
    for k in ks:
        for i in range(len(dataspecs_speclabel)):
            l[i] = l[i] + d[k][i]
    er_val_global = [_reduce_mean_parabola(it) for it in l]


    pulls = []
    for dataset_item in d:
        ditem = d[dataset_item]
        ds_pulls = []
        for tb, (alphas_global, error_global) in zip(ditem, er_val_global):
            if tb is None:
                ds_pulls.append(None)
                continue
            try:
                cv, error = _reduce_mean_parabola(tb)
            except ValueError:
                log.error(f"Concave parabola for {dataset_item}")
                ds_pulls.append(None)
            else:
                ds_pulls.append(_pulls_func(cv, alphas_global, error, error_global))

        pulls.append(ds_pulls)

    pulls = np.asarray(pulls,dtype=float)
    #Add total sum
    #sp = np.atleast_2d(np.nansum(pulls, axis=0))
    #pulls =np.concatenate([pulls, sp], axis=0)
    #fig,ax = barplot(pull_s.T, [*d.keys(), 'sum'], dataspecs_speclabel)

    fig,ax = barplot(pulls.T, d.keys(), dataspecs_speclabel)

    ax.legend()
    return fig

# Define aliases for functions with spelling mistakes in their names which have now been corrected
# Do this so that old runcards still work
plot_as_datasets_pseudorreplicas_chi2 = plot_as_datasets_pseudoreplicas_chi2
plot_dataspecs_pseudorreplica_means = plot_dataspecs_pseudoreplica_means

@figure
def plot_alphas_history(replica_alphaslog, number_alphas_history_to_plot="all"):
    number_of_replicas = len(replica_alphaslog)

    # Get history of the average alphas value
    max_epochs = max([len(i) for i in replica_alphaslog])
    arr = np.empty((number_of_replicas, max_epochs))
    for i in range(number_of_replicas):
        arr[i, :replica_alphaslog[i].size] = replica_alphaslog[i]
        arr[i,replica_alphaslog[i].size:] = replica_alphaslog[i][-1]
    # arr = np.ma.empty((number_of_replicas, max_epochs))
    # arr.mask = True
    # for i in range(number_of_replicas):
    #     arr[i, :replica_alphaslog[i].size] = replica_alphaslog[i]
    average_alphas_history = arr.mean(axis=0)

    # Get history of alphas for random replicas
    if number_alphas_history_to_plot == "all":
        number_alphas_history_to_plot = number_of_replicas 
    if number_of_replicas > number_alphas_history_to_plot:
        random_ordered_replica_indices = list(range(number_of_replicas))
        import random
        random.shuffle(random_ordered_replica_indices)
        random_replicas_to_plot = [random_ordered_replica_indices.pop() for _ in range(20)]
    else:
        random_replicas_to_plot = list(range(number_of_replicas))
    replica_alphaslog=[replica_alphaslog[i] for i in random_replicas_to_plot]

    # Make the figure
    fig, ax = plt.subplots()
    cmap = plt.cm.viridis
    aa=[i[0] for i in replica_alphaslog]
    norm = matplotlib.colors.Normalize(
        vmin=min(aa),
        vmax=max(aa)
        )
    for alphas_history in replica_alphaslog:
        ax.plot(alphas_history, color=cmap(norm(alphas_history[0])))
    ax.plot(average_alphas_history, color="red", label=r"mean $\alpha_s$ value")
    ax.set_xlabel("epochs")
    ax.set_ylabel(r"$\alpha_s$")
    ax.legend()
    if number_of_replicas > number_alphas_history_to_plot:
        ax.set_title(r"$\alpha_s$ history of "
                    f"{len(random_replicas_to_plot)} randomly drawn replica fits")
    else:
        ax.set_title(r"$\alpha_s$ history of the replica fits")
    return fig
