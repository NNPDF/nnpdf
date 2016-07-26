    # -*- coding: utf-8 -*-
"""
Utilities for reweighting studies.

Implements utilities for calculating the NNPDF weights and unweighted PDF sets.
It also allows for some basic statistics.
"""
import logging
from collections import OrderedDict


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import lhapdf
from reportengine.table import table
from reportengine.figure import figure
from reportengine.checks import make_check
from reportengine.formattingtools import spec_to_nice_name

from validphys.core import PDF
from validphys.results import abs_chi2_data, results
from validphys import checks
from validphys import lhaindex
from validphys.lhio import new_pdf_from_indexes

log = logging.getLogger(__name__)

__all__ = ('chi2_data_for_reweighting_experiments', 'make_unweighted_pdf',
           'nnpdf_weights', 'nnpdf_weights_numerator', 'p_alpha_study',
           'plot_p_alpha', 'reweighting_stats', 'unweighted_index',
           'make_pdf_from_filtered_outliers')

#TODO: implement this using reportengine expand when available
#use_t0 is to request that parameter to be set explicitly
@checks.check_pdf_is_montecarlo
def chi2_data_for_reweighting_experiments(reweighting_experiments, pdf, use_t0,
                                          t0set=None):
    """Like chi2data, but for reweighting experiments."""
    return [abs_chi2_data(results(exp,pdf,t0set,)) for exp in reweighting_experiments]


def nnpdf_weights_numerator(chi2_data_for_reweighting_experiments):
    """Compute the numerator of the NNPDF weights. This is useful for P(α),
    which uses a different normalization."""
    total_ndata = 0
    chi2s = np.zeros_like(chi2_data_for_reweighting_experiments[0][0].data)
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
    return np.exp(np.nansum(w*np.log(N/w))/N)

@table
def reweighting_stats(pdf, nnpdf_weights, p_alpha_study):
    """Compute varios statistics related to reweighting.

    Those are:
     - Number of initial replicas.
     - Effective number of replicas.
     - Median of the weightd.
     - The maximum value of P(alpha) in some sensible range.
    """
    er = effective_number_of_replicas(nnpdf_weights)
    initial_replicas = len(pdf) - 1
    median = np.median(nnpdf_weights)
    max_alpha = p_alpha_study.argmax()

    result = OrderedDict([
                          (r'N_{initial}', initial_replicas),
                          (r'$N_{eff}$', er),
                          (r'median($w$)', median),
                          (r'$max_{[0.5,4]}P(\alpha)$', max_alpha)
                         ])

    return pd.Series(result, index=result.keys())

def p_alpha_study(chi2_data_for_reweighting_experiments):
    """Compute P(α) in the range (0.5, 4)"""
    alphas = np.exp(np.linspace(np.log(0.5), np.log(4),31))
    vals = []
    for alpha in alphas:
        new_chi2 = [((type(res)(res.data/alpha)), central, ndata)
                    for (res,central,ndata) in
                    chi2_data_for_reweighting_experiments]
        new_ws = nnpdf_weights_numerator(new_chi2)
        val = np.sum(new_ws / alpha)
        vals.append(val)
    return pd.Series(np.array(vals), index=alphas)

@figure
def plot_p_alpha(p_alpha_study):
    """Plot the results of ``p_alpha_study``."""
    fig, ax = plt.subplots()
    ax.set_title(r"$P(\alpha)$")

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



@make_check
def _prepare_pdf_name(*, callspec, ns, environment, **kwargs):
    #TODO: Does this make any sense?
    if ns['output_path'] is not None:
        raise checks.CheckError("Output folder is not meant to be overwritten")
    output_path = environment.output_path / 'pdfsets'
    output_path.mkdir(exist_ok=True)
    ns['output_path'] = output_path

    set_name = ns['set_name']
    rootns = ns.maps[-1]
    if set_name is None:
        if 'nreplicas' in ns:
            suffix = ns['nreplicas']
        elif 'fit' in ns:
            suffix = ns['fit'].name
        else:
            raise RuntimeError("Bad namespace")
        set_name = spec_to_nice_name(rootns, callspec, str(suffix))
        ns['set_name'] = set_name

    if lhaindex.isinstalled(set_name):
        raise checks.CheckError("The PDF set that would be "
                         "generated already exists in the LHAPDF path:\n%s\n"
                         "Either delete it or explicitly assign a set_name for "
                         "the new PDF."
                         % lhaindex.finddir(set_name))

    #Ugly hack to allow analyzing the generated pdf some day (as in smpdf)
    if '_future_pdfs' not in rootns:
        rootns['_future_pdfs'] = {}

    future_pdfs = rootns['_future_pdfs']



    if set_name in future_pdfs:
        raise checks.CheckError("PDF set with name %s would already be "
                         "generated by another action and would be overwritten"
                         % set_name)

    lhapdf.pathsAppend(str(output_path))
    future_pdfs[set_name] = callspec


#TODO: This should return a PDF and nothing else. The environment should
#take care of saving the PDFs like we do for everything else.
@_prepare_pdf_name
@checks.check_can_save_grid
def make_unweighted_pdf(pdf, unweighted_index,
                        set_name:(str, type(None))=None, output_path=None,
                        installgrid:bool=True):
    """Generate an unweighted PDF set, from the prior ``pdf`` and the
    reweighting_experiments."""
    new_pdf_from_indexes(pdf=pdf, indexes=np.ravel(unweighted_index),
                         set_name=set_name, folder=output_path,
                         installgrid=installgrid)

    return PDF(set_name)

#Display this in the help
make_unweighted_pdf.highlight = 'pdfset'


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
@_prepare_pdf_name
def make_pdf_from_filtered_outliers(fit, replica_data, nsigma_cut:float,
                                    set_name:(str, type(None))=None,
                                    output_path=None,
                                    installgrid:bool=True):
    """Discard outliers with χ² to fitted data bigger than nsigma_cut and
    produce a new grid with the result."""
    chis = np.array([dt.chi2 for dt in replica_data])
    mean = np.mean(chis)
    std = np.std(chis)
    limit = mean + std*nsigma_cut
    #Replica indexes start at 1
    indexes = np.where(chis < limit)[0] + 1
    log.info("Mean ̉χ² is %.2f, and the threshold is %.2f. Discarded %d "
    "replicas out of %d.", mean, limit, len(chis) - len(indexes), len(chis))
    new_pdf_from_indexes(pdf=PDF(fit.name), indexes=indexes,
                         set_name=set_name, folder=output_path,
                         installgrid=installgrid)

make_pdf_from_filtered_outliers.highlight = 'pdfset'
