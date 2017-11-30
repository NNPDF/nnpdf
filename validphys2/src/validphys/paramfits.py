"""
paramfits.py

NOTE: This module is experimental and under develpment. The interfaces here
may change at will.

Functionality to determine parameters from a scan over
PDF fits. αs is so far the only example.

The functions here are high level and specialized, and rely on the more low
level modules (e.g. fitdata.py and results.py) for most of the functionality.

They also need to work around the limitations in libnnpdf, and so the
performance may not be optimal.
"""
import logging
import functools
from collections import namedtuple, defaultdict, Counter

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.mlab as mlab
import scipy.stats as stats

from reportengine.figure import figure, figuregen
from reportengine.table import table
from reportengine import collect
from reportengine.floatformatting import format_error_value_columns, format_number, ValueErrorTuple
from reportengine.checks import make_argcheck, CheckError, check_positive, check_not_empty
from NNPDF import pseudodata, single_replica, RandomGenerator

from validphys.core import PDF
from validphys.results import ThPredictionsResult, DataResult, chi2_breakdown_by_dataset
from validphys.plotutils import plot_horizontal_errorbars, marker_iter_plot, barplot, kde_plot

log = logging.getLogger(__name__)


PseudoReplicaExpChi2Data = namedtuple('PseudoReplicaChi2Data',
    ['experiment', 'dataset', 'ndata', 'chi2', 'nnfit_index'])

def computed_psedorreplicas_chi2(
        experiments, dataseed, pdf, fitted_replica_indexes,
        t0set:(PDF, type(None))):
    """Return a dataframe with the chi² of each replica wirh its corrsponding
    pseudodata (i.e. the one it was fitted with). The chi² is computed for
    both each experiment and each dataset in the experiment. The index of the
    dataframe is

    ['experiment', 'dataset', 'ndata' , 'nnfit_index']

    where 'experiment' is the name of the experiment, 'dataset' is the name of
    the dataset, or "Total" for the total value, 'ndata' is the corresponding
    number of points and 'nnfit_index' is the index specifying the
    pseudorreplica fluctuation.
    """

    #TODO: Everythning about this function is horrible. We need to rewrite
    #experiments.cc from scratch.

    #TODO: Do this somewhere else
    RandomGenerator.InitRNG(0,0)
    if t0set is not None:
        lt0 = t0set.load_t0()
    pdfname = pdf.name
    datas = []

    #No need to save these in the cache, so we call __wrapped__
    original_experiments = [e.load.__wrapped__(e) for e in experiments]
    sqrtcovmat_table = []
    log.debug("Generating dataset covmats")
    for exp in original_experiments:
        if t0set is not None:
            exp.SetT0(lt0)
        #The covariance matrices are currently very expensive to recompute.
        #Store them after computing T0
        sqrtcovmat_table.append([ds.get_sqrtcovmat() for ds in exp.DataSets()])

    for lhapdf_index, nnfit_index in enumerate(fitted_replica_indexes, 1):

        flutuated_experiments = pseudodata(original_experiments, dataseed, nnfit_index)
        lpdf = single_replica(pdfname, lhapdf_index)
        for expspec, exp, mats in zip(experiments, flutuated_experiments, sqrtcovmat_table):
            #We need to manage the memory
            exp.thisown = True

            th = ThPredictionsResult.from_convolution(pdf, expspec,
                loaded_data=exp, loaded_pdf=lpdf)


            results = DataResult(exp), th
            #The experiment already has T0. No need to set it again.
            #TODO: This is a hack. Get rid of this.
            chi2 = chi2_breakdown_by_dataset(results, exp, t0set=None,
                                             datasets_sqrtcovmat=mats)

            for i, (dsname,(value, ndata)) in enumerate(chi2.items()):
                data = PseudoReplicaExpChi2Data(
                    nnfit_index=nnfit_index,
                    experiment=expspec.name,
                    #We set the i so that the sort order is maintaned here.
                    dataset = (i, dsname),
                    ndata = ndata,
                    chi2=value
                    )
                datas.append(data)

    df =  pd.DataFrame(datas, columns=PseudoReplicaExpChi2Data._fields)
    df.set_index(['experiment', 'dataset', 'ndata' , 'nnfit_index'], inplace=True)
    df.sort_index(inplace=True)
    #Now that we have the order we like, we remove the i
    df.index.set_levels([x[1] for x in df.index.levels[1]], level=1, inplace=True)
    return df

#TODO: Probably fitcontext should set all of the variables required to compute
#this. But better setting
#them explicitly than setting some, se we require the user to do that.
fits_computed_psedorreplicas_chi2 = collect(computed_psedorreplicas_chi2, ('fits',))

dataspecs_computed_pseudorreplicas_chi2 = collect(computed_psedorreplicas_chi2, ('dataspecs',))

def _check_list_different(l, name):
    strs = [str(item) for item in l]
    if not len(set(strs))==len(l):
        counter = Counter(strs)
        duplicates = [k for k, v in counter.items() if v > 1]
        raise CheckError(f"{name} must be all different "
                         f"but there are duplicates: {duplicates}")

@make_argcheck
def _check_fits_different(fits):
    """Need this check because oterwise the pandas object gets confused"""
    return _check_list_different(fits, 'fits')

@make_argcheck
def _check_dataspecs_fits_different(dataspecs_fit):
    """Need this check because oterwise the pandas object gets confused"""
    return _check_list_different(dataspecs_fit, 'fits')

#TODO: Export the total here. Not having it is causing huge pain elsewhere.
@table
@_check_fits_different
def fits_matched_pseudorreplicas_chi2_table(fits, fits_computed_psedorreplicas_chi2):
    """Collect the chi^2 of the pseudoreplicas in the fits a single table,
    groped by nnfit_id.
    The columns come in two levels, fit name and (total chi², n).
    The indexes also come in two levels: nnfit_id and experiment name."""
    return pd.concat(fits_computed_psedorreplicas_chi2, axis=1, keys=map(str,fits))

@table
@_check_dataspecs_fits_different
def dataspecs_matched_pseudorreplicas_chi2_table(
        dataspecs_fit, dataspecs_computed_pseudorreplicas_chi2):
    """Like ``fits_matched_pseudorreplicas_chi2_table`` but for arbitrary dataspecs"""
    return fits_matched_pseudorreplicas_chi2_table(dataspecs_fit, dataspecs_computed_pseudorreplicas_chi2)


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

@make_argcheck
def _check_badcurves(badcurves):
    options = ['discard', 'minimum', 'allminimum']
    if badcurves not in  options:
        raise CheckError(f"badcurves must be one of {options}",
                         badcurves, options)




def _discard_sparse_curves(fits_replica_data_correlated,
        max_ndiscarded):
    """Return a table like  `fits_replica_data_correlated` where the replicas
    with too many discarded points have been filtered out."""

    df = fits_replica_data_correlated

    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded

    table = table[filt]
    return table, filt

@make_argcheck
def _check_discarded_string(max_ndiscarded):
    arg = max_ndiscarded
    if isinstance(arg,str):
        if arg != 'auto':
            raise CheckError("Expecting string to be 'auto'")

@_check_discarded_string
def discarded_mask(
    fits_replica_data_correlated_for_total,
    fits_as,
    max_ndiscarded:(int,str)='auto',
    autodiscard_confidence_level:float=0.99,
    trim_ndistant:int=0):

    """Return a table like  `fits_replica_data_correlated` where the replicas
    with too many discarded points have been filtered out.

    autodiscard_confidence_level is the student-T confidence level. Is normalised to 1
    and only is used if max_ndiscarded is set to 'auto'

    The automated discarding is done by estimating the uncertainty on the uncertainty by bootstrapping.

    The function returns a mask to be applied in fits_replica_data_with_discarded_replicas"""

    df = fits_replica_data_correlated_for_total[0]


    estimate = parabolic_as_determination(fits_as,df)
    best_as = np.mean(estimate)
    dist_best_as = -np.abs(best_as - fits_as)
    to_remove = np.argpartition(dist_best_as, trim_ndistant)[:trim_ndistant]
    as_mask = np.ones(df.shape[1], dtype=bool)
    as_mask[to_remove] = False



    if isinstance(max_ndiscarded,int):
        return _discard_sparse_curves(df,max_ndiscarded)[1], as_mask

    else:
        best_error = np.inf
        ndiscarded = range(len(fits_as),0,-1)
        for i in range(len(ndiscarded),0,-1):

            tablefilt_total, auto_filt = _discard_sparse_curves(df,ndiscarded[i-1])
            least_points = tablefilt_total.notnull().sum(axis=1).min()

            #Number of points that pass the cuts
            size = np.sum(auto_filt)

            #We can only fit a parabola with 3 points.
            #Use a fouth to have in principle some error estimate.
            if least_points > 3:
                parabolas = parabolic_as_determination(fits_as,tablefilt_total)
                bootstrap_est = np.random.choice(parabolas,(100000,size)).std(axis=1).std()
            else:
                bootstrap_est = np.inf

            stdT = stats.t.ppf((1-(1-autodiscard_confidence_level)/2), size-1)
            current_err = bootstrap_est*stdT

            if current_err < best_error:
                best_error = current_err
                best_filt = auto_filt

        return best_filt, as_mask

def fits_replica_data_with_discarded_replicas(
        discarded_mask,
        fits_replica_data_correlated):
    """Applies mask from discarded_mask to dataframes"""
    curve_mask, as_mask = discarded_mask

    discarded_replicas = fits_replica_data_correlated[curve_mask].copy()
    #Set these to Nan instead to masking them away in order to not break
    #all the apis that match this with fits_as.
    discarded_replicas.iloc[:,~as_mask] = np.NAN
    return discarded_replicas

def _get_parabola(asvals, chi2vals):
    chi2vals = np.ravel(chi2vals)
    filt =  np.isfinite(chi2vals)
    return np.polyfit(np.asarray(asvals)[filt], chi2vals[filt], 2)



def _parabolic_as_minimum_and_coefficient(fits_as,
        fits_replica_data_with_discarded_replicas,
        badcurves='discard'):
    """This is like parabolic_as_determination but without the as_transform
    functionality"""
    alphas = fits_as

    table = fits_replica_data_with_discarded_replicas.as_matrix()

    minimums = []
    quadratic = []
    asarr = np.asarray(alphas)
    for row in table:
        filt =  np.isfinite(row)
        if not filt.any():
            continue
        a,b,c = np.polyfit(asarr[filt], row[filt], 2)
        quadratic.append(a)
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
    quadratic = np.asarray(quadratic)
    minimums = np.asarray(minimums)
    return minimums, quadratic


def _parabolic_as_determination(fits_as,
        fits_replica_data_with_discarded_replicas,badcurves='discard'):

    return _parabolic_as_minimum_and_coefficient( fits_as,
                   fits_replica_data_with_discarded_replicas, badcurves)[0]

def quadratic_as_determination(fits_as,
        fits_replica_data_with_discarded_replicas,badcurves='discard'):

    return _parabolic_as_minimum_and_coefficient( fits_as,
                   fits_replica_data_with_discarded_replicas, badcurves)[1]



@make_argcheck
def _check_as_transform(as_transform):
    values = (None, 'log', 'exp')
    if not as_transform in values:
        raise CheckError(f"The allowed valued values for "
                         f"as_transform are {values}", str(as_transform),
                         values[1:])

@_check_badcurves
@_check_as_transform
def parabolic_as_determination(fits_as,
        fits_replica_data_with_discarded_replicas,
        badcurves='discard', as_transform:(str, type(None))=None):
    """Return the minima for alpha_s corresponding to the fitted curves.
    ``badcuves`` specifies what to do with concave replicas and can be one of
    'discard', 'allminimum'
    (which takes the minimum points
    for *all* the replicas without fitting a parabola) or
    'minimum' (which takes the minimum value for the concave replicas).

    as_transform can be None, 'log' or 'exp' and is applied to the as_values
    and then reversed for the minima.
    """
    if as_transform == 'log':
        fits_as = np.log(fits_as)
    elif as_transform == 'exp':
        fits_as = np.exp(fits_as)
    minimums = _parabolic_as_determination(
                   fits_as,
                   fits_replica_data_with_discarded_replicas, badcurves)
    if as_transform == 'log':
        minimums = np.exp(minimums)
    elif as_transform == 'exp':
        minimums = np.log(minimums)
    return minimums

def as_central_parabola(
        fits_as,
        fits_total_chi2):
    """Return the coefficients corresponding to the parabolic fit to the
    minimum of the pseudorreplicas"""
    return _get_parabola(fits_as, fits_total_chi2)

as_datasets_central_parabolas = collect(
        'as_central_parabola', ['fits_central_chi2_by_dataset_item'])

central_by_dataset_suptitle = collect('suptitle', ['fits_central_chi2_by_dataset_item'])
dataspecs_central_by_dataset_suptitle = collect('central_by_dataset_suptitle', ['dataspecs'])

central_by_dataset_ndata = collect(
    'ndata',
    ['fits_central_chi2_by_dataset_item',]
)
dataspecs_central_by_dataset_ndata = collect('central_by_dataset_ndata', ['dataspecs'])

by_dataset_suptitle = collect(
    'suptitle',
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)

dataspecs_dataset_suptitle = collect('by_dataset_suptitle', ['dataspecs'])


by_dataset_ndata = collect(
    'ndata',
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)

dataspecs_dataset_ndata = collect('by_dataset_ndata', ['dataspecs'])


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
    ax.plot(asarr, np.polyval(as_central_parabola, asarr)/ndata)

    best_as = np.mean(parabolic_as_determination_for_total)
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
    """Plot the cummulative total chi² for each of the datasets"""
    fig,ax  = plt.subplots()
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
                                     parabolic_as_determination_for_total):
    """Plot the cummulative difference between the χ² at the best global
    αs fit and the χ² at αs. If the difference is negative, it is set to zero.
    """
    fig,ax  = plt.subplots()
    nx = 100
    best_as = np.mean(parabolic_as_determination_for_total)
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
    ax.set_xlim(asarr[[0,-1]])

    return fig


@table
def derivative_dispersion_table(
                                as_datasets_central_parabolas,
                                fits_as,
                                central_by_dataset_suptitle,
                                as_determination_from_central_chi2_for_total):
    best_as = np.ravel(as_determination_from_central_chi2_for_total)[0]
    d = {}
    for label, p in zip(central_by_dataset_suptitle, as_datasets_central_parabolas):
        d[label] = np.polyval(np.polyder(p), best_as)

    s = pd.Series(d)
    s['SUM'] = np.sum(s)
    s['SUM QUADRATURE'] = np.sqrt(np.sum(s**2))
    res =  pd.DataFrame(s, columns=['Derivative'])

    return res

dataspecs_as_central_parabolas = collect('as_datasets_central_parabolas', ['dataspecs'])

def dataspecs_as_central_parabolas_map(
        dataspecs_speclabel,
        dataspecs_as_central_parabolas,
        dataspecs_central_by_dataset_suptitle,
        dataspecs_central_by_dataset_ndata):
    """Return a dict-like datastucture with the central chi² of the form:

        d[dataset_name][dataspec] = parabola_coefficients/ndata

    for all dataset items and dataspecs.
    """
    res = defaultdict(dict)
    for label, parabolas, dsnames, ndatas in zip(
            dataspecs_speclabel,
            dataspecs_as_central_parabolas,
            dataspecs_central_by_dataset_suptitle,
            dataspecs_central_by_dataset_ndata):
        for parabola, dsname, ndata  in zip(parabolas, dsnames, ndatas):
            res[dsname][label] = parabola/np.asarray(ndata)
    return res

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

def _aic(residuals, n, k):
    return 2*k + n*np.log(residuals) + 2*k*(k+1)/(n-k-1)

@table
def compare_aic(fits_as, fits_replica_data_with_discarded_replicas, suptitle):
    """Compare the  Akaike information criterion (AIC) for a
    parabolic and a cubic fit.  Note that
    this does **not** yield the actual AIC score, but only the piece
    neccessary to compare least squared fit (i.e. assuming
    iid gaussian noise for all points). This is:

        2*k + n*log(sum(residuals squared))

    The mean and standard deviation are taken across curves.
    Note that this always uses the *discard* criterion:
    That is, it ignores the curves that have no minimim."""
    alphas = fits_as
    asarr = np.asarray(alphas)

    aic2s = []
    aic3s = []

    table = fits_replica_data_with_discarded_replicas.as_matrix()
    for row in table:
        filt =  np.isfinite(row)
        asfilt = asarr[filt]
        rowfilt = row[filt]
        n = len(rowfilt)

        p2, res2, *stuff = np.polyfit(asfilt, rowfilt, 2, full=True)
        if p2[0] <= 0:
            pass
            #log.warning(f"Concave parabola computing AIC in {suptitle}")
        else:
            aic2 = _aic(res2, n, k=4)
            aic2s.append(aic2)

        p3, res3, *stuff = np.polyfit(asfilt, rowfilt, 3, full=True)

        extrema = np.roots(np.polyder(p3))
        #Cast away the zero complex part
        candidates = np.real(extrema[np.isreal(extrema)])
        if not len(candidates):
            pass
            #log.warning(f"Bad cubic minimum computing AIC in {suptitle}")
        else:
            aic3 = _aic(res3, n, k=5)
            aic3s.append(aic3)
    v2, e2 = np.mean(aic2s), np.std(aic2s)
    v3, e3 = np.mean(aic3s), np.std(aic3s)

    qp  = "Quadratic polynomial"
    cp = "Cubic polynomial"

    df = pd.DataFrame({'mean': {qp: v2, cp: v3}, 'error': {qp: e2, cp:e3},
                       'n minima':{qp: len(aic2s), cp: len(aic3s)}},
                      columns=['mean', 'error', 'n minima'])
    format_error_value_columns(df, 'mean', 'error', inplace=True)
    return df


def as_determination_from_central_chi2(fits_as, fits_total_chi2):
    """Return the alpha_s from the minimim chi² and the Delta_chi²=1 error
    from a quadratic fit to the total chi²."""
    alphas = fits_as
    chi2s = np.ravel(fits_total_chi2)
    a,b,c = np.polyfit(alphas, chi2s, 2)
    if a<=0:
        log.error("Found non convex parabola when computing the quadratic fit.")
        return np.nan, np.nan
    return ValueErrorTuple(-b/(2*a), 1/(np.sqrt(a)))

def parabolic_as_determination_with_tag(parabolic_as_determination, suptitle):
    """Convenience function to collect the arguments together. It is an identity"""
    return parabolic_as_determination, suptitle

def quadratic_as_determination_with_tag(quadratic_as_determination, suptitle):
    """Convenience function to collect the arguments together. It is an identity"""
    return quadratic_as_determination, suptitle

def as_determination_from_central_chi2_with_tag(
        as_determination_from_central_chi2, suptitle):
    """Convenience function to collect the arguments together. It is an identity"""

    return as_determination_from_central_chi2, suptitle


as_datasets_pseudorreplicas_chi2 = collect(
    parabolic_as_determination_with_tag,
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)

as_datasets_central_chi2 = collect(
    as_determination_from_central_chi2_with_tag,
    ['fits_central_chi2_by_dataset_item']
)

parabolic_as_determination_for_total = collect(parabolic_as_determination,
                                      ['matched_pseudorreplcias_for_total'])

as_determination_from_central_chi2_for_total = collect(
        as_determination_from_central_chi2, ['fits_central_chi2_for_total'])


quadratic_datasets_pseudorreplicas_chi2 = collect(
    quadratic_as_determination_with_tag,
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)


@figure
def plot_as_datasets_pseudorreplicas_chi2(as_datasets_pseudorreplicas_chi2):
    """Plot the error bars of the αs determination from pseudorreplicas
    by dataset item. Note that this only has meaning of preferred
    value for "Total", and the rest of the values are the minima of
    the partial χ²."""
    data, names = zip(*as_datasets_pseudorreplicas_chi2)
    cv, err = zip(*[(np.mean(dt), np.std(dt)) for dt in data])
    fig, ax = plot_horizontal_errorbars([cv], [err], names)
    ax.set_xlabel(r"$\alpha_S$")
    ax.set_title(r"$\alpha_S$ from pseudorreplicas")
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
def plot_as_datasets_compare(as_datasets_pseudorreplicas_chi2,
                             as_datasets_central_chi2,
                             marktotal:bool=True):
    """Plot the result of ``plot_as_datasets_pseudorreplicas_chi2`` and
    ``plot_as_exepriments_central_chi2`` together."""
    datapseudo, namespesudo = zip(*as_datasets_pseudorreplicas_chi2)
    cvpseudo, errpseudo = zip(*[(np.mean(dt), np.std(dt)) for dt in datapseudo])


    datacentral, namescentral = zip(*as_datasets_central_chi2)
    cvcentral, errcentral = zip(*datacentral)

    if namespesudo != namescentral:
        raise RuntimeError("Names do not coincide")

    fig, ax = plot_horizontal_errorbars(
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
def bootstrapping_stats_error(parabolic_as_determination, nresamplings:int=100000, suptitle=""):
    """Compute the bootstrapping uncertainty of the distribution of
    deterrminations of as, by resampling the list of points with replacement
    from the original sampling distribution `nresamplings` times
    and thn computing the standard deviation of the means."""
    distribution = parabolic_as_determination
    shape = (nresamplings, len(distribution))
    if not len(distribution):
        log.error("Cannot conpute stats error. Empty data.")
        return np.nan
    return np.random.choice(distribution, shape).mean(axis=1).std()


@check_positive('nresamplings')
def half_sample_stats_error(parabolic_as_determination, nresamplings:int=100000):
    """Like the bootstrapping error, but using only half og the data"""
    distribution = parabolic_as_determination[:len(parabolic_as_determination)//2]
    if not len(distribution):
        log.error("Cannot compute half stats. Too few data")
        return np.nan
    shape = (nresamplings, len(distribution))
    return np.random.choice(distribution, shape).mean(axis=1).std()




as_datasets_bootstrapping_stats_error = collect(bootstrapping_stats_error,
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)

as_datasets_half_sample_stats_error = collect(half_sample_stats_error,
    ['fits_matched_pseudorreplicas_chi2_by_dataset_item',]
)


#Don't write complicated column names everywhere
ps_mean = "pseudorreplica mean"
ps_error = "pseudorreplica error"
ps_stat_error = "pseudorreplica stat"
ps_half_stat_error = "pseudorreplica halfstat"
stats_ratio = r"$\frac{halfstat}{stat}/\sqrt 2$"
n = 'n'

stats_halfone = "cv selecting one half of the replicas"
err_halfone = "err selecting one half of the replicas"

stats_halfother = "cv selecting other half of the replicas"
err_halfonother = "err selecting other half of the replicas"

ps_cols = (ps_mean, ps_error ,n, ps_stat_error, ps_half_stat_error,
           stats_halfone, err_halfone, stats_halfother, err_halfonother)


cv_mean = "central mean"
cv_error = "central error"

def pseudorreplicas_stats_error(
        as_datasets_pseudorreplicas_chi2,
        as_datasets_bootstrapping_stats_error,
        as_datasets_half_sample_stats_error):
    """Return a dictionary (easily convertible to a DataFrame) with the mean,
    error and the measures of statistical error for each dataset."""
    d = defaultdict(dict)

    for (distribution, tag), statserr, halfstaterr in zip(
                as_datasets_pseudorreplicas_chi2,
                as_datasets_bootstrapping_stats_error,
                as_datasets_half_sample_stats_error):
        d[ps_mean][tag] = np.mean(distribution)
        d[n][tag] = len(distribution)
        d[ps_error][tag] = np.std(distribution)
        d[ps_stat_error][tag] = statserr
        d[ps_half_stat_error][tag] = halfstaterr
        d[stats_ratio][tag] = halfstaterr/statserr/np.sqrt(2)

        ldh = len(distribution)//2
        onehalf = distribution[:ldh]
        otherhalf = distribution[ldh:]
        d[stats_halfone][tag] = np.mean(onehalf)
        d[err_halfone][tag] = np.std(onehalf)
        d[stats_halfother][tag] = np.mean(otherhalf)
        d[err_halfonother][tag] = np.std(otherhalf)

    return dict(d)

datasepecs_pseudorreplica_stats_error = collect(pseudorreplicas_stats_error, ['dataspecs'])

#TODO: This is deprecated FAPP
@make_argcheck
def _check_dataset_items(dataset_items, dataspecs_dataset_suptitle):
    """Check that the dataset_items are legit."""
    if dataset_items is None:
        return
    try:
        s = set(dataset_items)
    except Exception as e:
        raise CheckError(f'dataset_items must be a list of strings: {e}') from e

    flat = [item for l in dataspecs_dataset_suptitle for item in l]
    d = s - set(flat)
    if d:
        raise CheckError(f"The floowing dataset_items are unrecognized: {d}")

@make_argcheck
def _check_speclabels_different(dataspecs_speclabel):
    """This is needed for grouping dataframes (and because
    generally indecated a bug)"""
    return _check_list_different(dataspecs_speclabel, 'dataspecs_speclabel')


def compare_determinations_table_impl(
        pseudorreplicas_stats_error,
        as_datasets_central_chi2):
    """Produce a table by experiment comparing the alpha_S determination
    from pseudorreplcias and from central values."""


    #Use this to get the right sorting
    d  = defaultdict(dict)
    tags = []
    for (cv, error), tag in as_datasets_central_chi2:
        d[cv_mean][tag] = cv
        d[cv_error][tag] = error

        tags.append(tag)

    d.update(pseudorreplicas_stats_error)

    df = pd.DataFrame(d, columns=[*ps_cols, cv_mean, cv_error])
    df = df.loc[tags]
    return df

@table
@_check_speclabels_different
def dataspecs_stats_error_table(
        datasepecs_pseudorreplica_stats_error,
        dataspecs_dataset_suptitle,
        dataspecs_speclabel,
        dataset_items:(type(None), list) = None,
        ):
    """Return a table with the stats errors of the pseudorreplica determination
    of each dataspec"""
    dfs = []
    for d in datasepecs_pseudorreplica_stats_error:
        df = pd.DataFrame(d, columns=ps_cols)
        format_error_value_columns(df, ps_mean, ps_error)
        format_error_value_columns(df, stats_halfone, err_halfone)
        format_error_value_columns(df, stats_halfother, err_halfonother)

        dfs.append(df)
    table =  pd.concat(dfs, axis=1, keys=dataspecs_speclabel)
    if dataset_items is not None:
        table = table.loc[dataset_items]
    return table

@table
def compare_determinations_table(compare_determinations_table_impl):
    """Return ``compare_determinations_table_impl`` formatted nicely"""
    df = compare_determinations_table_impl
    format_error_value_columns(df, ps_mean,
         ps_error, inplace=True)
    format_error_value_columns(df, cv_mean,
        cv_error, inplace=True)
    stats_cols = {ps_stat_error, ps_half_stat_error, stats_ratio}
    #Don't fail if/when we remove a table from here
    stats_cols &= set(df.columns)
    stats_cols = list(stats_cols)

    digits2 = functools.partial(format_number, digits=2)
    df[stats_cols] = df[stats_cols].applymap(digits2)
    return df

dataspecs_as_datasets_pseudorreplicas_chi2 = collect('as_datasets_pseudorreplicas_chi2', ['dataspecs'])

quad_as_datasets_pseudorreplicas_chi2 = collect('quadratic_datasets_pseudorreplicas_chi2',['dataspecs'])


#TODO: Deprecate fixup dataset_items earlier
@_check_speclabels_different
@_check_dataset_items
@table
def dataspecs_ndata_table(
            dataspecs_dataset_suptitle,
            dataspecs_dataset_ndata,
            dataspecs_speclabel,
            dataset_items:(list, type(None))=None):
    """Return a table with the same index as
    datasepecs_as_value_error_table_impl with the number of points
    per dataset."""
    d = {}
    for dslabel, datanames, ndatas in zip(dataspecs_speclabel,
                                          dataspecs_dataset_suptitle,
                                          dataspecs_dataset_ndata):
        d[dslabel] = dict(zip(datanames, ndatas))
    df = pd.DataFrame(d)
    if dataset_items is not None:
        df = df.loc[dataset_items]
    return df

@_check_speclabels_different
@_check_dataset_items
def datasepecs_quad_table_impl(
        quad_as_datasets_pseudorreplicas_chi2, dataspecs_speclabel,
        dataspecs_dataset_suptitle,
        dataset_items:(list, type(None)) = None,
        display_n:bool = False,
        ):
    """Return a table with the mean and error of the quadratic value of the parabolic
    determinations across dataspecs. If display_n is True, a column showing the number of points
    will be added to the table"""
    tables = []
    taglist = {}
    if display_n:
        cols = ['mean', 'error', 'n']
    else:
        cols = ['mean', 'error']
    for dets in quad_as_datasets_pseudorreplicas_chi2:
        d = defaultdict(dict)

        for distribution, tag in dets:
            d['mean'][tag] = np.mean(distribution)
            d['error'][tag] = np.std(distribution)

            if display_n:
                d['n'][tag] = len(distribution)
            taglist[tag] = None

        tables.append(pd.DataFrame(d, columns=cols))

    df = pd.concat(tables, axis=1, keys=dataspecs_speclabel)
    if dataset_items is None:
        ordered_keys = list(taglist)
    else:
        ordered_keys = dataset_items

    df = df.loc[ordered_keys]
    return df



@_check_speclabels_different
@_check_dataset_items
def datasepecs_as_value_error_table_impl(
        dataspecs_as_datasets_pseudorreplicas_chi2, dataspecs_speclabel,
        dataspecs_dataset_suptitle,
        dataset_items:(list, type(None)) = None,
        display_n:bool = False,
        ):
    """Return a table with the mean and error of the as determinations across
    dataspecs. If display_n is True, a column showing the number of points
    will be added to the table"""
    tables = []
    #Use the fact that in py3.6 a dict with None values is like an ordered set
    #TODO: A better way to build the dataframe?
    taglist = {}
    if display_n:
        cols = ['mean', 'error', 'n']
    else:
        cols = ['mean', 'error']
    for dets in dataspecs_as_datasets_pseudorreplicas_chi2:
        d = defaultdict(dict)

        for distribution, tag in dets:
            d['mean'][tag] = np.mean(distribution)
            d['error'][tag] = np.std(distribution)


            if display_n:
                d['n'][tag] = len(distribution)
            taglist[tag] = None

        tables.append(pd.DataFrame(d, columns=cols))


    df = pd.concat(tables, axis=1, keys=dataspecs_speclabel)
    if dataset_items is None:
        ordered_keys = list(taglist)
    else:
        ordered_keys = dataset_items

    df = df.loc[ordered_keys]



    return df

@table
def dataspecs_as_value_error_table(datasepecs_as_value_error_table_impl):
    """Return ``datasepecs_value_error_table_impl`` formatted nicely"""
    def f(x):
        return format_error_value_columns(x, x.columns[0], x.columns[1])
    return datasepecs_as_value_error_table_impl.groupby(level=0, axis=1).apply(f)

@table
def dataspecs_as_value_error_table_transposed(dataspecs_as_value_error_table):
    """Transposed version of ``dataspecs_as_value_error_table``.
    Useful for printing"""
    return dataspecs_as_value_error_table.T

@table
def dataspecs_quad_value_error_table(datasepecs_quad_table_impl):
    """Return ``datasepecs_value_error_table_impl`` formatted nicely"""
    def f(x):
        return format_error_value_columns(x, x.columns[0], x.columns[1])
    return datasepecs_quad_table_impl.groupby(level=0, axis=1).apply(f)



@figure
def plot_dataspecs_as_value_error(datasepecs_as_value_error_table_impl,
        dataspecs_fits_as,
        marktotal:bool=True, fix_limits:bool=True):
    """
    Plot the result for each dataspec of the pseudorreplica alpha_s
    determination based on the  partial chi² for each ``dataset_item``.

    If ``marktotal`` is True, a verical line will appear marking the position
    of the best fit.

    If ``fix_limits`` is True, the limits of the plot will span all the fitted
    values. Otherwise an heuristic will be used.

    """

    df = datasepecs_as_value_error_table_impl
    datalabels = df.columns.levels[0]
    catlabels = list(df.index)
    cvs = df.loc[:, (slice(None), 'mean')].T.as_matrix()
    errors = df.loc[:, (slice(None), 'error')].T.as_matrix()

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

@make_argcheck
def _check_first_is_total(fits_central_chi2_by_experiment_and_dataset):
    l = fits_central_chi2_by_experiment_and_dataset
    if not l or l[0]['experiment_label'] != 'Total':
        raise CheckError("Expecting that the first experiment is the total. You may need to set prepend_total=True")




@figure
@_check_first_is_total
def plot_as_value_error_central(as_datasets_central_chi2,
         marktotal:bool=True):
    """Plot the result of ``plot_as_datasets_pseudorreplicas_chi2`` and
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
    ax.plot(x, kde_pulls(x), label="Kernal Density Estimation of pulls")
    ax.hist(pulls,normed=True,bins=4)
    ax.grid(False)
    ax.plot(x, mlab.normpdf(x, mean_pulls, std_dev),label="Normalised gaussian fit")
    ax.legend()

    return fig

@figure
def plot_pull_plots_global_min(datasepecs_as_value_error_table_impl,
        dataspecs_fits_as,dataspecs_speclabel,hide_total:bool=True):

    """Plots the pulls of individual experiments as a barplot."""

    df = datasepecs_as_value_error_table_impl
    tots_error = df.loc['Total', (slice(None), 'error')].as_matrix()
    tots_mean = df.loc['Total', (slice(None), 'mean')].as_matrix()

    if hide_total:
        df = df.loc[df.index != 'Total']

    catlabels = list(df.index)
    cvs = df.loc[:, (slice(None), 'mean')].as_matrix()
    errors = df.loc[:, (slice(None), 'error')].as_matrix()

    pulls = _pulls_func(cvs,tots_mean,errors,tots_error).T

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
    datasepecs_as_value_error_table_impl,
    datasepecs_quad_table_impl,
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
    df = datasepecs_as_value_error_table_impl
    df2 = datasepecs_quad_table_impl


    tots_mean = df.loc['Total', (slice(None), 'mean')].as_matrix()

    if hide_total:
        df = df.loc[df.index != 'Total']
        df1 = df1.loc[df1.index != 'Total']
        df2 = df2.loc[df2.index != 'Total']


    cvs = df.loc[:, (slice(None), 'mean')].T.as_matrix()
    quad_weights = df2.loc[:, (slice(None), 'mean')].T.as_matrix()

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
def plot_pull_gaussian_fit_pseudo(datasepecs_as_value_error_table_impl,
        dataspecs_fits_as,dataspecs_speclabel,hide_total:bool=True):

    """Bins the pulls computed in pull_plots_global_min and overlays
    the normalised gaussian fit and KDE to the histogram of pulls"""

    df = datasepecs_as_value_error_table_impl
    tots_error = df.loc['Total', (slice(None), 'error')].T.as_matrix()
    tots_mean = df.loc['Total', (slice(None), 'mean')].T.as_matrix()

    if hide_total:
        df = df.loc[df.index != 'Total']

    cvs = df.loc[:, (slice(None), 'mean')].T.as_matrix()
    errors = df.loc[:, (slice(None), 'error')].T.as_matrix()

    for label, i in zip(dataspecs_speclabel, range(len(cvs))):
        pulls = _pulls_func(cvs[i],tots_mean[i],errors[i],tots_error[i])

        mean_pulls = np.mean(pulls)
        std_dev = np.std(pulls)
        x = np.linspace(min(pulls),max(pulls), 100)

        kde_pulls = stats.gaussian_kde(pulls, bw_method='silverman')
        fig, ax = plt.subplots()

        #ax.set_title(f"Histogram of pulls for {label} dataset")
        ax.set_xlabel(r"Pull")
        ax.plot(x, kde_pulls(x), label="Kernal Density Estimation of pulls")
        ax.hist(pulls,normed=True,bins=4)
        ax.grid(False)
        ax.plot(x, mlab.normpdf(x, mean_pulls, std_dev),label="Normalised gaussian fit")
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

dataspecs_fits_as = collect('fits_as', ['dataspecs'])

by_dataset_as_chi2 = collect(
        fits_replica_data_with_discarded_replicas,
        ['fits_matched_pseudorreplicas_chi2_by_dataset_item',])

dataspecs_fits_replica_data_with_discarded_replicas = collect(
        'by_dataset_as_chi2', ['dataspecs'])

@check_not_empty('dataspecs_dataset_suptitle')
def dataspecs_chi2_by_dataset_dict(dataspecs_dataset_suptitle,
                        dataspecs_fits_replica_data_with_discarded_replicas,
                        dataspecs_fits_as,
                        ):
    """Return a table-like dict with the
    suptitle: [<list of tables>]

    where each table is ``fits_replica_data_with_discarded_replicas`` resolved
    for the given dataset in each of the dataspecs.
    """
    allkeys = set(dataspecs_dataset_suptitle[0])
    for newkeys in dataspecs_dataset_suptitle[1:]:
        newkeys = set(newkeys)
        symdiff = newkeys ^ allkeys
        if symdiff:
            log.warning(f"Some datasets are not "
                        f"common across all dataspecs {symdiff}")
            allkeys |= symdiff

    res = defaultdict(list)
    for keys, values, asvals in zip(dataspecs_dataset_suptitle,
            dataspecs_fits_replica_data_with_discarded_replicas,
            dataspecs_fits_as):
        for k, v in zip(keys, values):
            v.columns = asvals
            res[k].append(v)
        for k in allkeys-set(keys):
            res[k].append(None)
    return res

@figuregen
@_check_dataset_items
def plot_dataspecs_pseudorreplica_means(
        dataspecs_chi2_by_dataset_dict,
        dataspecs_speclabel,
        dataset_items:(list, type(None))=None,
        ):
    """ Plot the mean chi² from data to pseudorreplica, over replicas in a fit
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
@_check_dataset_items
@check_positive('examples_per_item')
def plot_dataspecs_parabola_examples(
        dataspecs_chi2_by_dataset_dict,
        dataspecs_speclabel,
        dataset_items:(list, type(None))=None,
        examples_per_item:int = 2,
        random_seed:int = 0,
        ):
    """Sample ``examples_per_item`` replica_indexes for each of the
    ``dataset_items``. Yield a plot with the parabolic fit, as resoved for
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
                y = vals.as_matrix()
                ax.plot(asvals, y, **next(im), label=label,
                         color=color, linestyle='none', lw=0.5)
                a,b,c = parabola = _get_parabola(asvals, y)
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

    distribution = parabolic_as_determination

    fig, ax = plt.subplots()

    kde_plot(distribution)
    ax.legend()
    ax.set_title(f"{suptitle}")
    ax.set_xlabel(r"$\alpha_S$")
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

@table
def export_fits_computed_psedorreplicas_chi2(fits_computed_psedorreplicas_chi2):
    """Hack to force writting the CSV output"""
    return fits_computed_psedorreplicas_chi2
