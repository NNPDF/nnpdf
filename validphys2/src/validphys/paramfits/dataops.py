"""
dataops.py

This module implements the core functionality of the paramfits package,
currently focused on the αs determination.
It computes various statistics and formats the data in a form suitable to
be consumed by plotting functions.
"""
import logging
import warnings
import functools
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.stats as stats

from reportengine import collect
from reportengine.floatformatting import format_error_value_columns, ValueErrorTuple, format_number
from reportengine.checks import make_argcheck, CheckError, check_positive, check_not_empty
from reportengine.table import table

from validphys.checks import (
    check_fits_different,
    check_dataspecs_fits_different,
    check_speclabels_different,
)

log = logging.getLogger(__name__)

class StandardSampleWrapper:
    """A class that holds a flat array of data, and has 'location' and
    'scale' properties, which are the mean and the standard deviation.
    The purpose of the class is to explicitly disallow calling np.mean
    and np.std on the result while preserving functional backward
    compatibility with the runcards."""
    def __init__(self, data):
        self.data = data

    @property
    def location(self):
        return np.mean(self.data)

    @property
    def scale(self):
        return np.std(self.data)

class RobustSampleWrapper:
    """Similar to ``StandardSampleWrapper``, but location and scale
    are implemented as the median and the 68% interval respectively."""
    def __init__(self, data):
        self.data = data

    @property
    def location(self):
        return np.median(self.data)

    @property
    def scale(self):
        return np.asscalar(np.diff(np.percentile(self.data, [15.87, 84.13])))/2


def get_parabola(asvals, chi2vals):
    """Return the three coefficients of a parabola χ²(αs) given a set of
    asvals and a set of χ² values.
    Only the finite χ² values are taken into account."""
    chi2vals = np.ravel(chi2vals)
    filt =  np.isfinite(chi2vals)
    return np.polyfit(np.asarray(asvals)[filt], chi2vals[filt], 2)

#TODO: Export the total here. Not having it is causing huge pain elsewhere.
@table
@check_fits_different
def fits_matched_pseudoreplicas_chi2_table(fits, fits_computed_pseudoreplicas_chi2):
    """Collect the chi^2 of the pseudoreplicas in the fits a single table,
    groped by nnfit_id.
    The columns come in two levels, fit name and (total chi², n).
    The indexes also come in two levels: nnfit_id and experiment name."""
    return pd.concat(fits_computed_pseudoreplicas_chi2, axis=1, keys=map(str,fits))

@table
@check_dataspecs_fits_different
def dataspecs_matched_pseudoreplicas_chi2_table(
        dataspecs_fit, dataspecs_computed_pseudoreplicas_chi2):
    """Like ``fits_matched_pseudoreplicas_chi2_table`` but for arbitrary dataspecs"""
    return fits_matched_pseudoreplicas_chi2_table(dataspecs_fit, dataspecs_computed_pseudoreplicas_chi2)

@make_argcheck
def _check_badcurves(badcurves):
    options = ['discard', 'minimum', 'allminimum']
    if badcurves not in options:
        raise CheckError(f"badcurves must be one of {options}",
                         badcurves, options)




def _discard_sparse_curves(fits_replica_data_correlated,
        max_ndiscarded):
    """Return a table like `fits_replica_data_correlated` where the replicas
    with too many discarded points have been filtered out."""

    df = fits_replica_data_correlated

    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1, min_count=1) < max_ndiscarded

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

    """Return a table like `fits_replica_data_correlated` where the replicas
    with too many discarded points have been filtered out.

    autodiscard_confidence_level is the student-T confidence level. Is normalised to 1
    and only is used if max_ndiscarded is set to 'auto'

    The automated discarding is done by estimating the uncertainty on the uncertainty by bootstrapping.

    The function returns a mask to be applied in fits_replica_data_with_discarded_replicas"""

    df = fits_replica_data_correlated_for_total[0]

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        estimate = parabolic_as_determination(fits_as,df)
    best_as = estimate.location
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
            least_points = tablefilt_total.notnull().sum(axis=1, min_count=1).min()

            #Number of points that pass the cuts
            size = np.sum(auto_filt)

            #We can only fit a parabola with 3 points.
            #Use a fouth to have in principle some error estimate.
            if least_points > 3:
                #Apply as mask before fitting
                tablefilt_total = tablefilt_total.copy()
                tablefilt_total.iloc[:,~as_mask] = np.NAN
                parabolas = parabolic_as_determination(fits_as,tablefilt_total).data
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


def _parabolic_as_minimum_and_coefficient(fits_as,
        fits_replica_data_with_discarded_replicas,
        badcurves='discard'):
    """This implements the fitting in ``parabolic_as_determination``.
    Returns the minimum and the set of locations."""
    alphas = fits_as

    table = fits_replica_data_with_discarded_replicas.values

    minimums = []
    quadratic = []
    asarr = np.asarray(alphas)
    for row in table:
        filt = np.isfinite(row)
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

    return _parabolic_as_minimum_and_coefficient(fits_as,
                   fits_replica_data_with_discarded_replicas, badcurves)[1]



@make_argcheck
def _check_as_transform(as_transform):
    values = (None, 'log', 'exp', 'logshift')
    if not as_transform in values:
        raise CheckError(f"The allowed valued values for "
                         f"as_transform are {values}", str(as_transform),
                         values[1:])

@make_argcheck
def _check_parabolic_as_statistics(parabolic_as_statistics):
    values = 'standard', 'robust'
    if parabolic_as_statistics not in values:
        raise CheckError('The allowed values for'
                f'`parabolic_as_statistics` are {values}',
                parabolic_as_statistics, values)


@_check_badcurves
@_check_as_transform
@_check_parabolic_as_statistics
def parabolic_as_determination(
        fits_as,
        fits_replica_data_with_discarded_replicas,
        badcurves='discard',
        as_transform:(str, type(None))=None,
        parabolic_as_statistics:str='standard'):
    """Return the minima for alpha_s corresponding to the fitted curves.
    ``badcuves`` specifies what to do with concave replicas and can be one of
    'discard', 'allminimum'
    (which takes the minimum points
    for *all* the replicas without fitting a parabola) or
    'minimum' (which takes the minimum value for the concave replicas).

    If ``parabolic_as_statistics`` is ``"standard"``, means and standard
    deviations will be used to compute statstics. Otherwise, if it is
    ``"robust"``, medians and 68% intervals will be used.

    as_transform can be None, 'log', 'logshift' (``log(1+αs)``) or
    'exp' and is applied to the as_values and then reversed for the
    minima.
    """
    if as_transform == 'log':
        fits_as = np.log(fits_as)
    elif as_transform == 'logshift':
        fits_as = np.log(fits_as + 1)
    elif as_transform == 'exp':
        fits_as = np.exp(fits_as)
    minimums = _parabolic_as_determination(
                   fits_as,
                   fits_replica_data_with_discarded_replicas, badcurves)
    # Invert the transform
    if as_transform == 'log':
        minimums = np.exp(minimums)
    elif as_transform == 'exp':
        minimums = np.log(minimums)
    elif as_transform == 'logshift':
        minimums = np.exp(minimums - 1)
    if parabolic_as_statistics == 'standard':
        res = StandardSampleWrapper(minimums)
    elif parabolic_as_statistics == 'robust':
        res = RobustSampleWrapper(minimums)
    else:
        raise RuntimeError("Unknown `parabolic_as_statistics`.")
    return res

def as_central_parabola(
        fits_as,
        fits_total_chi2):
    """Return the coefficients corresponding to the parabolic fit to the
    minimum of the pseudoreplicas"""
    return get_parabola(fits_as, fits_total_chi2)

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
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

dataspecs_dataset_suptitle = collect('by_dataset_suptitle', ['dataspecs'])


by_dataset_ndata = collect(
    'ndata',
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

dataspecs_dataset_ndata = collect('by_dataset_ndata', ['dataspecs'])


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
    res = pd.DataFrame(s, columns=['Derivative'])

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
        for parabola, dsname, ndata in zip(parabolas, dsnames, ndatas):
            res[dsname][label] = parabola/np.asarray(ndata)
    return res

def _aic(residuals, n, k):
    return 2*k + n*np.log(residuals) + 2*k*(k+1)/(n-k-1)

@table
def compare_aic(fits_as, fits_replica_data_with_discarded_replicas, suptitle):
    """Compare the Akaike information criterion (AIC) for a
    parabolic and a cubic fit. Note that
    this does **not** yield the actual AIC score, but only the piece
    necessary to compare least squared fit (i.e. assuming
    iid gaussian noise for all points). This is:

        2*k + n*log(sum(residuals squared))

    The mean and standard deviation are taken across curves.
    Note that this always uses the *discard* criterion:
    That is, it ignores the curves that have no minimum."""
    alphas = fits_as
    asarr = np.asarray(alphas)

    aic2s = []
    aic3s = []

    table = fits_replica_data_with_discarded_replicas.values
    for row in table:
        filt = np.isfinite(row)
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

    qp = "Quadratic polynomial"
    cp = "Cubic polynomial"

    df = pd.DataFrame({'mean': {qp: v2, cp: v3}, 'error': {qp: e2, cp:e3},
                       'n minima':{qp: len(aic2s), cp: len(aic3s)}},
                      columns=['mean', 'error', 'n minima'])
    format_error_value_columns(df, 'mean', 'error', inplace=True)
    return df


def as_determination_from_central_chi2(fits_as, fits_total_chi2):
    """Return the alpha_s from the minimum chi² and the Delta_chi²=1 error
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


as_datasets_pseudoreplicas_chi2 = collect(
    parabolic_as_determination_with_tag,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

as_datasets_central_chi2 = collect(
    as_determination_from_central_chi2_with_tag,
    ['fits_central_chi2_by_dataset_item']
)

parabolic_as_determination_for_total = collect(parabolic_as_determination,
                                      ['matched_pseudoreplicas_for_total'])

as_determination_from_central_chi2_for_total = collect(
        as_determination_from_central_chi2, ['fits_central_chi2_for_total'])


quadratic_datasets_pseudoreplicas_chi2 = collect(
    quadratic_as_determination_with_tag,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

dataspecs_parabolic_as_determination_for_total = collect(
        'parabolic_as_determination_for_total', ['dataspecs'])

@check_positive('nresamplings')
def bootstrapping_stats_error(parabolic_as_determination, nresamplings:int=100000, suptitle=""):
    """Compute the bootstrapping uncertainty of the distribution of
    determinations of as, by resampling the list of points with replacement
    from the original sampling distribution `nresamplings` times
    and then computing the standard deviation of the means."""
    distribution = parabolic_as_determination.data
    shape = (nresamplings, len(distribution))
    if not len(distribution):
        log.error("Cannot compute stats error. Empty data.")
        return np.nan
    return np.random.choice(distribution, shape).mean(axis=1).std()

@check_positive('nresamplings')
def bootstrapping_stats_error_on_the_error(parabolic_as_determination, nresamplings:int=100000, suptitle=""):
    """Compute the bootstrapping uncertainty of standard deviation on the parabolic determination."""
    distribution = parabolic_as_determination.data
    shape = (nresamplings, len(distribution))
    if not len(distribution):
        log.error("Cannot compute stats error. Empty data.")
        return np.nan
    return np.random.choice(distribution, shape).std(axis=1).std()


@check_positive('nresamplings')
def half_sample_stats_error(parabolic_as_determination, nresamplings:int=100000):
    """Like the bootstrapping error, but using only half og the data"""
    sample = parabolic_as_determination.data
    distribution = sample[:len(sample)//2]
    if not len(distribution):
        log.error("Cannot compute half stats. Too few data")
        return np.nan
    shape = (nresamplings, len(distribution))
    return np.random.choice(distribution, shape).mean(axis=1).std()




as_datasets_bootstrapping_stats_error = collect(bootstrapping_stats_error,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

as_datasets_bootstrapping_stats_error_on_the_error = collect(
    bootstrapping_stats_error_on_the_error,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

as_datasets_half_sample_stats_error = collect(half_sample_stats_error,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)


#Don't write complicated column names everywhere
ps_mean = "pseudoreplica mean"
ps_error = "pseudoreplica error"
ps_stat_error = "pseudoreplica stat"
ps_half_stat_error = "pseudoreplica halfstat"
stats_ratio = r"$\frac{halfstat}{stat}/\sqrt 2$"
n = 'n'

stats_halfone = "cv selecting one half of the replicas"
err_halfone = "err selecting one half of the replicas"

stats_halfother = "cv selecting other half of the replicas"
err_halfonother = "err selecting other half of the replicas"

stats_err_err = "Stat error on the error"

ps_cols = (ps_mean, ps_error ,n, ps_stat_error, ps_half_stat_error,
           stats_halfone, err_halfone, stats_halfother, err_halfonother,
           stats_err_err)


cv_mean = "central mean"
cv_error = "central error"

def pseudoreplicas_stats_error(
        as_datasets_pseudoreplicas_chi2,
        as_datasets_bootstrapping_stats_error,
        as_datasets_bootstrapping_stats_error_on_the_error,
        as_datasets_half_sample_stats_error):
    """Return a dictionary (easily convertible to a DataFrame) with the mean,
    error and the measures of statistical error for each dataset."""
    d = defaultdict(dict)

    for (distribution, tag), statserr, staterrerr, halfstaterr in zip(
                as_datasets_pseudoreplicas_chi2,
                as_datasets_bootstrapping_stats_error,
                as_datasets_bootstrapping_stats_error_on_the_error,
                as_datasets_half_sample_stats_error):
        d[ps_mean][tag] = distribution.location
        d[n][tag] = len(distribution.data)
        d[ps_error][tag] = distribution.scale
        d[ps_stat_error][tag] = statserr
        d[ps_half_stat_error][tag] = halfstaterr
        d[stats_ratio][tag] = halfstaterr/statserr/np.sqrt(2)
        d[stats_err_err][tag] = staterrerr

        ldh = len(distribution.data)//2
        onehalf = distribution.data[:ldh]
        otherhalf = distribution.data[ldh:]
        d[stats_halfone][tag] = np.mean(onehalf)
        d[err_halfone][tag] = np.std(onehalf)
        d[stats_halfother][tag] = np.mean(otherhalf)
        d[err_halfonother][tag] = np.std(otherhalf)

    return dict(d)

dataspecs_pseudoreplica_stats_error = collect(pseudoreplicas_stats_error, ['dataspecs'])

#TODO: This is deprecated FAPP
@make_argcheck
def check_dataset_items(dataset_items, dataspecs_dataset_suptitle):
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
        raise CheckError(f"The following dataset_items are unrecognized: {d}")



def compare_determinations_table_impl(
        pseudoreplicas_stats_error,
        as_datasets_central_chi2):
    """Produce a table by experiment comparing the alpha_S determination
    from pseudoreplicas and from central values."""


    #Use this to get the right sorting
    d  = defaultdict(dict)
    tags = []
    for (cv, error), tag in as_datasets_central_chi2:
        d[cv_mean][tag] = cv
        d[cv_error][tag] = error

        tags.append(tag)

    d.update(pseudoreplicas_stats_error)

    df = pd.DataFrame(d, columns=[*ps_cols, cv_mean, cv_error])
    df = df.loc[tags]
    return df

@table
@check_speclabels_different
def dataspecs_stats_error_table(
        dataspecs_pseudoreplica_stats_error,
        dataspecs_dataset_suptitle,
        dataspecs_speclabel,
        dataset_items:(type(None), list) = None,
        ):
    """Return a table with the stats errors of the pseudoreplica determination
    of each dataspec"""
    dfs = []
    for d in dataspecs_pseudoreplica_stats_error:
        df = pd.DataFrame(d, columns=ps_cols)
        format_error_value_columns(df, ps_mean, ps_error)
        format_error_value_columns(df, stats_halfone, err_halfone)
        format_error_value_columns(df, stats_halfother, err_halfonother)

        dfs.append(df)
    table = pd.concat(dfs, axis=1, keys=dataspecs_speclabel)
    if dataset_items is not None:
        table = table.loc[dataset_items]
    return table

@table
def dataspecs_stats_error_table_transposed(dataspecs_stats_error_table):
    """Transposed version of dataspecs_stats_error_table for display
    purposes."""
    return dataspecs_stats_error_table.T

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

dataspecs_as_datasets_pseudoreplicas_chi2 = collect('as_datasets_pseudoreplicas_chi2', ['dataspecs'])

quad_as_datasets_pseudoreplicas_chi2 = collect('quadratic_datasets_pseudoreplicas_chi2',['dataspecs'])


#TODO: Deprecate fixup dataset_items earlier
@check_speclabels_different
@check_dataset_items
@table
def dataspecs_ndata_table(
            dataspecs_dataset_suptitle,
            dataspecs_dataset_ndata,
            dataspecs_speclabel,
            dataset_items:(list, type(None))=None):
    """Return a table with the same index as
    dataspecs_as_value_error_table_impl with the number of points
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

@check_speclabels_different
@check_dataset_items
def dataspecs_quad_table_impl(
        quad_as_datasets_pseudoreplicas_chi2, dataspecs_speclabel,
        dataspecs_dataset_suptitle,
        dataset_items:(list, type(None)) = None,
        display_n:bool = False,
        ):
    """Return a table with the mean and error of the quadratic coefficient of the parabolic
    determinations across dataspecs. If display_n is True, a column showing the number of points
    will be added to the table"""
    tables = []
    taglist = {}
    if display_n:
        cols = ['mean', 'error', 'n']
    else:
        cols = ['mean', 'error']
    for dets in quad_as_datasets_pseudoreplicas_chi2:
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



@check_speclabels_different
@check_dataset_items
def dataspecs_as_value_error_table_impl(
        dataspecs_as_datasets_pseudoreplicas_chi2, dataspecs_speclabel,
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
    for dets in dataspecs_as_datasets_pseudoreplicas_chi2:
        d = defaultdict(dict)

        for distribution, tag in dets:
            d['mean'][tag] = distribution.location
            d['error'][tag] = distribution.scale


            if display_n:
                d['n'][tag] = len(distribution.data)
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
def dataspecs_as_value_error_table(dataspecs_as_value_error_table_impl):
    """Return ``dataspecs_value_error_table_impl`` formatted nicely"""
    def f(x):
        return format_error_value_columns(x, x.columns[0], x.columns[1])
    return dataspecs_as_value_error_table_impl.groupby(level=0, axis=1).apply(f)

@table
def dataspecs_as_value_error_table_transposed(dataspecs_as_value_error_table):
    """Transposed version of ``dataspecs_as_value_error_table``.
    Useful for printing"""
    return dataspecs_as_value_error_table.T

@table
def dataspecs_quad_value_error_table(dataspecs_quad_table_impl):
    """Return ``dataspecs_value_error_table_impl`` formatted nicely"""
    def f(x):
        return format_error_value_columns(x, x.columns[0], x.columns[1])
    return dataspecs_quad_table_impl.groupby(level=0, axis=1).apply(f)

dataspecs_fits_as = collect('fits_as', ['dataspecs'])

by_dataset_as_chi2 = collect(
        fits_replica_data_with_discarded_replicas,
        ['fits_matched_pseudoreplicas_chi2_by_dataset_item',])

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


as_dataset_pseudodata = collect(
    fits_replica_data_with_discarded_replicas,
    ['fits_matched_pseudoreplicas_chi2_by_dataset_item',]
)

@table
def as_parabolic_coefficient_table(
        fits_as,
        by_dataset_suptitle,
        as_dataset_pseudodata):
    """Return a table of the parabolic fit of each dataset item, for each
    correlated replica. The index is the correlated_replica index and there
    are four columns for each dataset: 'a', 'b' and 'c' corresponding to the
    parabolic coefficients and 'min', which is ``-b/2/a`` if 'a' is positive,
    and NaN otherwise."""
    alphas = np.asarray(fits_as)
    tb_polys = []
    #Easier to repeat the polyfit code here than to change the format of the
    #various outputs.
    for tb in as_dataset_pseudodata:
        polys = []
        for row in np.asarray(tb):
            filt =  np.isfinite(row)
            if not filt.any():
                polys.append(np.array([np.nan]*3))
                continue
            poly = np.polyfit(alphas[filt], row[filt], 2)
            polys.append(poly)
        frame = pd.DataFrame(polys, index=tb.index, columns=['a','b', 'c'])
        frame['min'] = -frame['b']/2/frame['a']
        frame.loc[frame['a']<0,'min'] = np.nan
        tb_polys.append(frame)
    final_table = pd.concat(tb_polys, axis=1, keys=by_dataset_suptitle)
    return final_table

# Define aliases for functions with spelling mistakes in their names which have now been corrected
# Do this so that old runcards still work
fits_matched_pseudorreplicas_chi2_table = fits_matched_pseudoreplicas_chi2_table
dataspecs_matched_pseudorreplicas_chi2_table = dataspecs_matched_pseudoreplicas_chi2_table
as_datasets_pseudorreplicas_chi2 = as_datasets_pseudoreplicas_chi2
quadratic_datasets_pseudorreplicas_chi2 = quadratic_datasets_pseudoreplicas_chi2
pseudorreplicas_stats_error = pseudoreplicas_stats_error
datasepecs_pseudorreplica_stats_error = dataspecs_pseudoreplica_stats_error
dataspecs_as_datasets_pseudorreplicas_chi2 = dataspecs_as_datasets_pseudoreplicas_chi2
quad_as_datasets_pseudorreplicas_chi2 = quad_as_datasets_pseudoreplicas_chi2
datasepecs_quad_table_impl = dataspecs_quad_table_impl
datasepecs_as_value_error_table_impl = dataspecs_as_value_error_table_impl
