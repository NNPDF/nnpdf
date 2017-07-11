"""
paramfits.py

Functionality to determine parameters from a scan over
PDF fits. αs is so far the only example.

The functions here are high level and specialized, and rely on the more low
level modules (e.g. fitdata.py and results.py) for most of the functionality.

They also need to work around the limitations in libnnpdf, and so the
performance may not be optimal.
"""
import logging
import itertools
from collections import namedtuple, defaultdict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from reportengine.figure import figure
from reportengine.table import table
from reportengine import collect
from reportengine.checks import make_argcheck, CheckError, check_positive
from NNPDF import pseudodata, single_replica, RandomGenerator

from validphys.core import PDF
from validphys.results import ThPredictionsResult, DataResult, abs_chi2_data_experiment
from validphys.plotutils import expand_margin, marker_iter_plot, centered_range

log = logging.getLogger(__name__)


PseudoReplicaExpChi2Data = namedtuple('PseudoReplicaChi2Data',
    ['nnfit_index', 'experiment', 'central_chi2', 'ndata'])

def computed_psedorreplicas_chi2(experiments, dataseed, pdf,
                                 fitted_replica_indexes, t0set:PDF):
    """Return the chi2 the chi² of the pseudodata"""

    #TODO: Everythning about this function is horrible. We need to rewrite
    #experiments.cc from scratch.

    #TODO: Do this somewhere else
    RandomGenerator.InitRNG(0,0)

    lt0 = t0set.load_t0()
    pdfname = pdf.name
    datas = []

    original_experiments = [e.load.__wrapped__(e) for e in experiments]
    for exp in original_experiments:
        exp.SetT0(lt0)

    for lhapdf_index, nnfit_index in enumerate(fitted_replica_indexes, 1):
        #No need to save these in the cache, so we call __wrapped__

        flutuated_experiments = pseudodata(original_experiments, dataseed, nnfit_index)
        lpdf = single_replica(pdfname, lhapdf_index)
        for expspec, exp in zip(experiments, flutuated_experiments):
            #We need to manage the memory
            exp.thisown = True

            th = ThPredictionsResult.from_convolution(pdf, expspec,
                loaded_data=exp, loaded_pdf=lpdf)


            results = DataResult(exp), th
            chi2 = abs_chi2_data_experiment(results)

            data = PseudoReplicaExpChi2Data(
                nnfit_index=nnfit_index,
                experiment=expspec.name,
                central_chi2=chi2.central_result,
                ndata=chi2.ndata)
            datas.append(data)
    df =  pd.DataFrame(datas, columns=PseudoReplicaExpChi2Data._fields)
    df.set_index(['nnfit_index', 'experiment'], inplace=True)
    return df

#TODO: Probably fitcontext should set all of the variables required to compute
#this. But better setting
#them explicitly than setting some, se we require the user to do that.
fits_computed_psedorreplicas_chi2 = collect(computed_psedorreplicas_chi2, ('fits',))

@make_argcheck
def _check_fits_different(fits):
    """Need this check because oterwise the pandas object gets confused"""
    if not len(set(map(str, fits)))==len(fits):
        raise CheckError("fits must be all different but there are duplicates.")

@table
@_check_fits_different
def fits_matched_pseudorreplicas_chi2_table(fits, fits_computed_psedorreplicas_chi2):
    """Collect the chi^2 of the pseudoreplicas in the fits a single table,
    groped by nnfit_id.
    The columns come in two levels, fit name and (total chi², n).
    The indexes also come in two levels: nnfit_id and experiment name."""
    return pd.concat(fits_computed_psedorreplicas_chi2, axis=1, keys=map(str,fits))




@figure
def plot_fits_as_profile(fits_pdfs, fits_total_chi2, suptitle=None):
    """Plot the total central chi² as a function of the value of α_s.
    Note that this plots as a function of the key "AlphaS_MZ" in the LHAPDF
    file, which is annoyingly *not* α_s(MZ) for Nf<5."""
    fig, ax = plt.subplots()
    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    #Could be a transposed data frame
    fits_total_chi2 = np.ravel(fits_total_chi2)
    ax.plot(alphas, fits_total_chi2)
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²/N_{dat}$')
    if suptitle:
        fig.suptitle(suptitle)
    return fig

@make_argcheck
def _check_badcurves(badcurves):
    options = ['discard', 'minimum', 'allminimum']
    if badcurves not in  options:
        raise CheckError(f"badcurves must be one of {options}",
                         badcurves, options)


def fits_replica_data_with_discarded_replicas(fits_replica_data_correlated,
        max_ndiscarded:int=4):
    df = fits_replica_data_correlated
    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded
    table = table[filt]
    #table = np.atleast_2d(np.mean(table, axis=0))
    return table


@_check_badcurves
def parabolic_as_determination(fits_pdfs,
        fits_replica_data_with_discarded_replicas,
        badcurves='discard'):
    """Return the minima for alpha_s corresponding to the fitted curves."""
    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]

    table = fits_replica_data_with_discarded_replicas.as_matrix()

    minimums = []
    asarr = np.asarray(alphas)
    for row in table:
        filt =  np.isfinite(row)
        a,b,c = np.polyfit(asarr[filt], row[filt], 2)
        if badcurves == 'allminimum':
            minimums.append(asarr[filt][np.argmin(row[filt])])
        elif a>0:
            minimums.append(-b/2/a)
        elif badcurves == 'discard':
            pass
        elif badcurves == 'minimum':
            minimums.append(asarr[filt][np.argmin(row[filt])])
        else:
            raise RuntimeError("Unknown bad curves.")

    minimums = np.asarray(minimums)
    return minimums


def as_determination_from_central_chi2(fits_pdfs, fits_total_chi2):
    """Return the alpha_s from the minimim chi² and the Delta_chi²=1 error
    from a quadratic fit to the total chi²."""
    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    chi2s = np.ravel(fits_total_chi2)
    a,b,c = np.polyfit(alphas, chi2s, 2)
    if a<=0:
        log.error("Found non convex parabola when computing the quadratic fit.")
        return np.nan, np.nan
    return -b/(2*a), 1/(np.sqrt(a))


def parabolic_as_determination_with_tag(parabolic_as_determination, suptitle):
    """Convenience function to collect the arguments together. It is an identity"""
    return parabolic_as_determination, suptitle


def as_determination_from_central_chi2_with_tag(
        as_determination_from_central_chi2, suptitle):
    """Convenience function to collect the arguments together. It is an identity"""

    return as_determination_from_central_chi2, suptitle

as_exepriments_psudoreplicas_chi2 = collect(
    parabolic_as_determination_with_tag,
    ['fits_matched_pseudorreplicas_chi2_output_by_experiment']
)

as_experiments_central_chi2 = collect(
    as_determination_from_central_chi2_with_tag,
    ['fits_absolute_chi2_output_by_experiment',]
)

#TODO: Move to plotutils
def _plot_horizontal_error_bars(cvs, errors, categorylabels, datalabels=None):
    w,h = plt.rcParams["figure.figsize"]
    rescale = max(1, 1 + 0.1*(len(categorylabels) - 7))
    fig, ax = plt.subplots(figsize=(w, h*rescale))
    if datalabels is None:
        datalabels = itertools.repeat(None)
    y = np.arange(len(categorylabels))
    ax.yaxis.set_ticks(y)
    ax.yaxis.set_ticklabels(categorylabels)
    mi = marker_iter_plot()


    distance = 0.5/len(cvs)
    pos = centered_range(len(cvs), distance=distance)

    for cv, err, lb, markerspec, shift in zip(cvs, errors, datalabels, mi, pos):
        ax.errorbar(cv, y+shift, xerr=err, linestyle='none', label=lb,
                    **markerspec)
    ax.set_xlim(*expand_margin(np.percentile(cvs, 15),
                               np.percentile(cvs, 85),
                               1.1))
    return fig, ax


@figure
def plot_as_exepriments_psudoreplicas_chi2(as_exepriments_psudoreplicas_chi2):
    """Plot the error bas of the alha_s determination from pseudorreplicas
    by experiment"""
    data, names = zip(*as_exepriments_psudoreplicas_chi2)
    cv, err = zip(*[(np.mean(dt), np.std(dt)) for dt in data])
    fig, ax = _plot_horizontal_error_bars([cv], [err], names)
    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ from pseudorreplicas")
    return fig

@figure
def plot_as_exepriments_central_chi2(as_experiments_central_chi2):
    """Plot the error bas of the alha_s determination from centrla chi²
    by experiment"""
    data, names = zip(*as_experiments_central_chi2)
    cv, err = zip(*data)
    fig, ax = _plot_horizontal_error_bars([cv], [err], names)
    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ from central chi²")
    return fig


#TODO: take compare_determinations_table
@figure
def plot_as_exepriments_compare(as_exepriments_psudoreplicas_chi2, as_experiments_central_chi2, marktotal:bool=True):
    """Plot the result of ``plot_as_exepriments_psudoreplicas_chi2`` and
    ``plot_as_exepriments_central_chi2`` together."""
    datapseudo, namespesudo = zip(*as_exepriments_psudoreplicas_chi2)
    cvpseudo, errpseudo = zip(*[(np.mean(dt), np.std(dt)) for dt in datapseudo])


    datacentral, namescentral = zip(*as_experiments_central_chi2)
    cvcentral, errcentral = zip(*datacentral)

    if namespesudo != namescentral:
        raise RuntimeError("Names do not coincide")

    fig, ax = _plot_horizontal_error_bars(
        [cvcentral, cvpseudo], [errcentral, errpseudo], namescentral,
        [r'Central $\chi^2$', r'Pseudorreplica $\chi^2$']
    )
    if marktotal:
        try:
            pos = namespesudo.index('Total')
        except ValueError:
            log.error("Asked to mark total, but it was not provided.")
        else:
            ax.axvline(cvcentral[pos], color='C0', linewidth=0.5, linestyle='--')
            ax.axvline(cvpseudo[pos], color='C1', linewidth=0.5, linestyle='--')

    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ determination")
    ax.legend()
    return fig

@check_positive('nresamplings')
def bootstrapping_stats_error(parabolic_as_determination, nresamplings:int=100000):
    """Compute the bootstrapping uncertainty of the distribution of
    deterrminations of as, by resampling the list of points with replacement
    from the original sampling distribution `nresamplings` times
    and thn computing the standard deviation of the means."""
    distribution = parabolic_as_determination
    shape = (nresamplings, len(distribution))
    return np.random.choice(distribution, shape).mean(axis=1).std()

@check_positive('nresamplings')
def half_sample_stats_error(parabolic_as_determination, nresamplings:int=100000):
    """Like the bootstrapping error, but using only half og the data"""
    distribution = parabolic_as_determination[:len(parabolic_as_determination)//2]
    shape = (nresamplings, len(distribution))
    return np.random.choice(distribution, shape).mean(axis=1).std()


as_experiments_bootstrapping_stats_error = collect(bootstrapping_stats_error,
    ['fits_matched_pseudorreplicas_chi2_output_by_experiment'],
)

as_experiments_half_sample_stats_error = collect(half_sample_stats_error,
    ['fits_matched_pseudorreplicas_chi2_output_by_experiment'],
)



@table
def compare_determinations_table(as_exepriments_psudoreplicas_chi2,
                                 as_experiments_central_chi2,
                                 as_experiments_bootstrapping_stats_error,
                                 as_experiments_half_sample_stats_error):
    """Produce a table by experiment comparing the alpha_S determination
    from pseudorreplcias and from central values."""
    d = defaultdict(dict)

    ps_mean = "pseudirreplica mean"
    ps_error = "pseudorreplica error"
    ps_stat_error = "pseudorreplica stat"
    ps_half_stat_error = "pseudorreplica halfstat"
    stats_ratio = r"$\frac{halfstat}{stat}/\sqrt 2$"

    cv_mean = "central mean"
    cv_error = "centeal error"

    for (distribution, tag), statserr, halfstaterr in zip(
                as_exepriments_psudoreplicas_chi2,
                as_experiments_bootstrapping_stats_error,
                as_experiments_half_sample_stats_error):
        d[ps_mean][tag] = np.mean(distribution)
        d[ps_error][tag] = np.std(distribution)
        d[ps_stat_error][tag] = statserr
        d[ps_half_stat_error][tag] = halfstaterr
        d[stats_ratio][tag] = halfstaterr/statserr/np.sqrt(2)


    for (cv, error), tag in as_experiments_central_chi2:
        d[cv_mean][tag] = cv
        d[cv_error][tag] = error

    return pd.DataFrame(d, columns=[ps_mean, ps_error,
        ps_stat_error, ps_half_stat_error, stats_ratio,
        cv_mean, cv_error])



@figure
def plot_fitted_replicas_as_profiles_matched(fits_pdfs,
        fits_replica_data_with_discarded_replicas,
        parabolic_as_determination, suptitle=None):
    """Plot chi²(as) keeping the replica nnfit index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """
    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]


    minimums = parabolic_as_determination

    table = fits_replica_data_with_discarded_replicas.as_matrix()

    fig, ax = plt.subplots()

    from matplotlib.collections import LineCollection

    lc = LineCollection([list(x for x in zip(alphas, t) if np.isfinite(x[1])) for t in table])
    lc.set_array(minimums)
    lc.set_clim(*np.percentile(minimums, (5,95)))
    ax.add_collection(lc)
    ax.set_xlim(min(alphas), max(alphas))
    ax.set_ylim(np.nanmin(table), np.nanmax(table))
    fig.colorbar(lc, label=r"Preferred $\alpha_S$")
    ax.set_title(rf"$\alpha_S$ = ${np.mean(minimums):.4f} \pm {np.std(minimums):.4f}$ N={len(minimums)}")
    if suptitle:
        fig.suptitle(suptitle)

    #ax.plot(alphas, np.array(table).T, color='#ddddcc')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi^2$')
    return fig


@figure
def plot_poly_as_fit(fits_pdfs,
        fits_replica_data_correlated, max_ndiscarded:int=4, polorder:int=2,
        suptitle=None):
    """Plot a polynomial fit of chi²(as) of `degree polorder`, keeping the
    replica index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """

    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    df = fits_replica_data_correlated
    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded
    table = table[filt]
    table = table.as_matrix()
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


