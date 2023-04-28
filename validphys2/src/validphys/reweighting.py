    # -*- coding: utf-8 -*-
"""
Utilities for reweighting studies.

Implements utilities for calculating the NNPDF weights and unweighted PDF sets.
It also allows for some basic statistics.
"""
import logging
import warnings
from collections import OrderedDict


import numpy as np
import pandas as pd
from scipy import optimize

from reportengine import collect
from reportengine.table import table
from reportengine.figure import figure

from validphys.core import PDF, Filter
from validphys import checks
from validphys.lhio import new_pdf_from_indexes
from validphys.pdfoutput import pdfset
from validphys.dataplots import plot_training_validation
from validphys import plotutils

log = logging.getLogger(__name__)

__all__ = ('chi2_data_for_reweighting_experiments', 'make_unweighted_pdf',
           'nnpdf_weights', 'nnpdf_weights_numerator', 'p_alpha_study',
           'plot_p_alpha', 'reweighting_stats', 'unweighted_index',
           'make_pdf_from_filtered_outliers')


chi2_data_for_reweighting_experiments_inner = collect(
    'abs_chi2_data_experiment', ['reweighting_experiments']
)

# This is to add checks and to request that use_t0 is set explicitly
# pylint: disable=unused-argument
@checks.check_pdf_is_montecarlo
def chi2_data_for_reweighting_experiments(
    chi2_data_for_reweighting_experiments_inner, use_t0
):
    return chi2_data_for_reweighting_experiments_inner


def nnpdf_weights_numerator(chi2_data_for_reweighting_experiments):
    """Compute the numerator of the NNPDF weights. This is useful for P(α),
    which uses a different normalization."""
    total_ndata = 0
    chi2s = np.zeros_like(chi2_data_for_reweighting_experiments[0][0].error_members())
    for data in chi2_data_for_reweighting_experiments:
        res, _, ndata = data
        total_ndata += ndata
        chi2s += res.data

    chi2s = np.ravel(chi2s)

    logw = ((total_ndata - 1)/2)*np.log(chi2s) - 0.5*chi2s
    logw -= np.max(logw)
    w = np.exp(logw)
    return w

@table
#will call list[0]
@checks.check_not_empty('reweighting_experiments')
def nnpdf_weights(chi2_data_for_reweighting_experiments):
    """Compute the replica weights according to the NNPDF formula."""
    numerator = nnpdf_weights_numerator(chi2_data_for_reweighting_experiments)
    return pd.DataFrame(numerator/np.sum(numerator),
                        index=np.arange(1, len(numerator) + 1))

def effective_number_of_replicas(w):
    N = len(w)
    w = w*N/np.sum(w)
    return np.exp(np.nansum(w*np.log(N)-w*np.log(w))/N)

@table
def reweighting_stats(pdf, nnpdf_weights, p_alpha_study):
    """Compute various statistics related to reweighting.

    Those are:
     - Number of initial replicas.
     - Effective number of replicas.
     - Median of the weightd.
     - The maximum value of P(alpha) in some sensible range.
    """
    er = effective_number_of_replicas(nnpdf_weights)
    initial_replicas = pdf.get_members()
    median = np.median(nnpdf_weights)
    max_alpha = p_alpha_study.idxmax()

    result = OrderedDict([
                          (r'$N_{\mathrm{initial}}$', initial_replicas),
                          (r'$N_{eff}$', er),
                          (r'median($w$)', median),
                          (r'$max_{[0.5,4]}P(\alpha)$', max_alpha)
                         ])

    return pd.Series(result, index=result.keys(), name='Reweighting stats')

def _get_p_alpha_val(alpha, chi2_data_for_reweighting_experiments):
     new_chi2 = [((type(res)(res.error_members()/alpha**2)), central, ndata)
                    for (res,central,ndata) in
                    chi2_data_for_reweighting_experiments]

     new_ws = nnpdf_weights_numerator(new_chi2)
     val = np.sum(new_ws / alpha)
     return val


def _get_p_alpha_vals(alphas, chi2_data_for_reweighting_experiments):
    vals = []
    for alpha in alphas:
        val = _get_p_alpha_val(alpha, chi2_data_for_reweighting_experiments)
        vals.append(val)

    return vals


def p_alpha_study(chi2_data_for_reweighting_experiments):
    """Compute P(α) in an automatic range"""



    small = 0.5

    f = lambda alpha: -_get_p_alpha_val(alpha, chi2_data_for_reweighting_experiments)
    #Ignore warnings for nonsensical alphas
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        alphamax=  optimize.fmin(f, 10, disp=False)[0]

    """
    #First find the maximum
    small = 0.5

    limclose=4
    #No, 100 is not a crazy big value
    limfar = 100

    nclose = 80
    nfar = 40

    alphas = np.zeros(nclose+nfar)
    alphas[:nclose] = np.linspace(small, limclose, nclose, endpoint=False)
    alphas[nclose:] = np.linspace(limclose,limfar,nfar)
    vals = _get_p_alpha_vals(alphas, chi2_data_for_reweighting_experiments)
    big = 1.2*alphas[np.argmax(vals)]

    #Compute more points around the maximum
    alphas = np.linspace(small, big, 100)
    """

    big = 1.4*alphamax


    alphas = np.linspace(small, big, 1000)
    vals = _get_p_alpha_vals(alphas, chi2_data_for_reweighting_experiments)


    return pd.Series(np.array(vals), index=alphas)

@figure
def plot_p_alpha(p_alpha_study):
    """Plot the results of ``p_alpha_study``."""
    fig, ax = plotutils.subplots()
    ax.set_title(r"$P(\alpha)$")

    xmax = p_alpha_study.idxmax()
    ymax = p_alpha_study[xmax]
    ax.axvline(xmax, color='red', linestyle='--')
    ax.annotate(r'$\alpha=%.2f$'%xmax, (xmax,(ymax-ax.get_ylim()[0])/2), )

    ax.set_yticklabels([])
    ax.set_xlabel(r'$\alpha$')


    ax.plot(p_alpha_study)
    return fig

@table
def unweighted_index(nnpdf_weights, nreplicas:int=100):
    """The index of the input replicas that corresponds to an unweighted set,
    for the given weights. This can be saved for testing purposes."""
    nnpdf_weights = np.ravel(nnpdf_weights)
    res = 1 + np.random.choice(len(nnpdf_weights), size=nreplicas, p=nnpdf_weights)
    return pd.DataFrame(res, index=np.arange(1,nreplicas+1))

@pdfset
@checks.check_can_save_grid
def make_unweighted_pdf(pdf, unweighted_index,
                        set_name:str, output_path=None,
                        installgrid:bool=True):
    """Generate an unweighted PDF set, from the prior ``pdf`` and the
    reweighting_experiments. The PDF is written to a `pdfsets` directory of
    the output folder. Return the relative path of the newly created PDF."""


    new_pdf_from_indexes(pdf=pdf, indexes=np.ravel(unweighted_index),
                         set_name=set_name, folder=output_path,
                         installgrid=installgrid)

@checks.make_check
def _check_cut(ns, *args, **kwargs):
    cut = ns['nsigma_cut']
    msg = "'nsigma_cut' must be a float greater than zero."
    try:
        if cut>0:
            return
    #TODO: Check types for parameters automatically
    except TypeError as e:
        raise checks.CheckError(msg) from e
    raise checks.CheckError(msg)

@checks.check_has_fitted_replicas
@_check_cut
def chi2filtered_index(fit, replica_data, nsigma_cut:float):
    """Discard outliers until the χ² to fitted data of all replicas has
       a smaller deviation than `nsigma_cut`
       (in units of the standard deviation of the final ensemble)
    """

    chis = np.array([dt.chi2 for dt in replica_data])
    newchis = chis
    oldmean = np.mean(newchis)
    while True:
        mean = np.mean(newchis)
        std = np.std(newchis)
        limit = mean + std*nsigma_cut
        #Replica indexes start at 1
        indexes = np.where(chis < limit)[0]


        if len(newchis) == len(indexes):
            break


        newchis = chis[indexes]

    log.info("%s: Mean ̉χ² was %.2f, and the threshold is %.2f. Discarded %d "
    "replicas out of %d. Now the mean is %.2f.", fit.label, oldmean, limit,
    len(chis) - len(indexes), len(chis), np.mean(newchis))

    label = '$\chi^2$ Filter: %.2f'%nsigma_cut


    return Filter(indexes, label, nsigma_cut=nsigma_cut)

def negative_filtered_index(count_negative_points, filter_Q=75):
    cut = np.percentile(count_negative_points, filter_Q)
    indexes =  np.where(count_negative_points <= cut)[0]
    label = 'Positivity Selection Q=%.2f' % filter_Q
    log.info("Positivity cut for Q=%.2f is at %.0f negative points" % (filter_Q, cut))
    return Filter(indexes, label, filter_Q=filter_Q)

@pdfset
@checks.check_can_save_grid
def make_pdf_from_filtered_outliers(fit, chi2filtered_index,
                                    set_name:str,
                                    output_path=None,
                                    installgrid:bool=True):
    """Produce a new grid with the result of chi2filtered_index"""

    indexes = chi2filtered_index.indexes + 1 #libnnpdf nonsense
    new_pdf_from_indexes(pdf=PDF(fit.name), indexes=indexes,
                         set_name=set_name, folder=output_path,
                         installgrid=installgrid)

#TODO: Use the filter framework here when it exists
@figure
def plot_chi2filtered_training_validation(fit, nsigma_cut, replica_data, chi2filtered_index):
    """Like `plot_training_validation`, but apply `chi2filtered_index` mask."""
    return plot_training_validation(fit, replica_data, dict([chi2filtered_index.as_pair()]))

#TODO: Use the filter framework here when it exists
@figure
def plot_posfiltered_training_validation(fit, replica_data, negative_filtered_index):
    """Like `plot_training_validation`, but apply `chi2filtered_index` mask."""
    return plot_training_validation(fit, replica_data, dict([negative_filtered_index.as_pair()]))


p_alpha_for_all_datasets = collect(p_alpha_study, ('reweight_all_datasets',))

chi_2_for_all_datasets = collect(chi2_data_for_reweighting_experiments,
                                 ('reweight_all_datasets',))

@table
def p_alpha_all_datasets_table(p_alpha_for_all_datasets, reweight_all_datasets, chi_2_for_all_datasets):
    """Compute and display P(alpha) and chi² for all datasets
    in all experiments."""
    data = []
    chilabel = r'$\chi^2$'
    modelabel = r'mode P($\alpha$)'
    for series, rexp, chis in zip(p_alpha_for_all_datasets, reweight_all_datasets, chi_2_for_all_datasets):
        central = chis[0][1]/chis[0][2]
        exp = rexp['reweighting_experiments'][0]
        data.append({'name':exp.datasets[0].name,
                     modelabel:series.argmax(),
                     chilabel: central,
                    })
    return pd.DataFrame.from_records(data, index='name',
                                     columns=('name', modelabel, chilabel))
