# -*- coding: utf-8 -*-
"""
Tools for computing and plotting effective exponents.
"""
from __future__ import generator_stop

import logging
import warnings
import numbers
import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.figure import figuregen
from reportengine.table import table
from reportengine.floatformatting import format_number, significant_digits
from reportengine.compat import yaml

from validphys.checks import check_positive, check_pdf_normalize_to, make_argcheck, check_xlimits
from validphys.pdfplots import BandPDFPlotter, PDFPlotter
from validphys.pdfbases import check_basis, Basis
from validphys.core import PDF, FitSpec


import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)

@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def alpha_eff(pdf: PDF, *,
              xmin: numbers.Real = 1e-6,
              xmax: numbers.Real = 1e-3,
              npoints: int = 200,
              Q: numbers.Real = 1.65,
              basis: (str, Basis),
              flavours: (list, tuple, type(None)) = None):
    """Return a list of xplotting_grids containing the value of the effective
    exponent alpha at the specified values of x and flavour.
    alpha is relevant at small x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    #Loading the filter map of the fit/PDF
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    if npoints == 2:
        xGrid = np.array([xmin, xmax])
    else:
        xGrid = pdfgrids.xgrid(xmin, xmax, 'log', npoints)

    pdfGrid = pdfgrids.xplotting_grid(
        pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
    pdfGrid_values = pdfGrid.grid_values
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        alphaGrid_values = -np.log(abs(pdfGrid_values/xGrid))/np.log(xGrid)
        alphaGrid_values[alphaGrid_values == - np.inf] = np.nan  # when PDF_i =0
    alphaGrid = pdfGrid._replace(grid_values=alphaGrid_values)
    return alphaGrid

@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def beta_eff(pdf, *,
             xmin: numbers.Real = 0.6,
             xmax: numbers.Real = 0.9,
             npoints: int = 200,
             Q: numbers.Real = 1.65,
             basis: (str, Basis),
             flavours: (list, tuple, type(None)) = None):
    """Return a list of xplotting_grids containing the value of the effective
    exponent beta at the specified values of x and flavour.
    beta is relevant at large x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    if npoints == 2:
        xGrid = np.array([xmin, xmax])
    else:
        xGrid = pdfgrids.xgrid(xmin, xmax, 'linear', npoints)


    pdfGrid = pdfgrids.xplotting_grid(
        pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
    pdfGrid_values = pdfGrid.grid_values
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        betaGrid_values = np.log(abs(pdfGrid_values/xGrid))/np.log(1-xGrid)
        betaGrid_values[betaGrid_values == -np.inf] = np.nan  # when PDF_i =0
    betaGrid = pdfGrid._replace(grid_values=betaGrid_values)

    return betaGrid  # .grid_values

class PreprocessingPlotter(PDFPlotter):
    """ Class inherenting from BandPDFPlotter, changing title and ylabel to reflect the effective
    exponent being plotted.
    """

    def __init__(self, exponent, *args,  **kwargs):
        self.exponent = exponent
        super().__init__(*args, **kwargs)

    def get_title(self, parton_name):
        return fr"$\{self.exponent}_e$ for ${parton_name}$ at {format_number(self.Q, 3)} Gev"

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return "Ratio to {}".format(self.normalize_pdf.label)
        else:
            return fr"$\{self.exponent}_e$ for ${parton_name}$"

def get_alpha_lines(effective_exponents_table):
    """Given an effective_exponents_table returns a dictionary which for `prev` and `next` contains
    the bounds of the alpha effective exponent for all flavours, used to plot horizontal lines on
    the alpha effective exponent plots.
    """
    df = effective_exponents_table[0]
    limits = {}
    limits['prev'] = df.iloc[0::2, [0, 1]].values
    limits['next'] = df.iloc[0::2, [2, 3]].values
    return limits

def get_beta_lines(effective_exponents_table):
    """Same as `get_alpha_lines` but for beta"""
    df = effective_exponents_table[0]
    limits = {}
    limits['prev'] = df.iloc[1::2, [0, 1]].values
    limits['next'] = df.iloc[1::2, [2, 3]].values
    return limits

fits_alpha_lines = collect('get_alpha_lines', ('fits',))
fits_beta_lines = collect('get_beta_lines', ('fits',))

class ExponentBandPlotter(BandPDFPlotter, PreprocessingPlotter):
    def __init__(self, hlines, exponent, fits, *args,  **kwargs):
        super().__init__(exponent, *args,  **kwargs)
        self.hlines = hlines
        self.fits = fits

    def draw(self, pdf, grid, flstate):
        # Here we assume each pdf has a corresponding fit, which is true by construction.
        errdown, errup = super().draw(pdf, grid, flstate)
        pdf_index = self.pdfs.index(pdf)
        hlines = self.hlines[pdf_index]
        xmin, xmax = flstate.ax.get_xlim()
        handle = flstate.ax.hlines(
            hlines['prev'][flstate.flindex, :], xmin=xmin, xmax=xmax, linestyle='-.')
        label = f'prev ({self.fits[pdf_index].label})'
        flstate.handles.append(handle)
        flstate.labels.append(label)

        handle = flstate.ax.hlines(
            hlines['next'][flstate.flindex, :], xmin=xmin, xmax=xmax, linestyle=':')
        label = f'next ({self.fits[pdf_index].label})'
        flstate.handles.append(handle)
        flstate.labels.append(label)
        # need to return xgrid shaped object but with hlines taken into account to get plots nice
        new_errdown = min(
            [*errdown, *hlines['next'][flstate.flindex, :], *hlines['prev'][flstate.flindex, :]])
        new_errup = max(
            [*errup, *hlines['next'][flstate.flindex, :], *hlines['prev'][flstate.flindex, :]])
        return new_errdown*np.ones_like(errdown), new_errup*np.ones_like(errup)


alpha_eff_pdfs = collect('alpha_eff', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_alphaEff_internal(
        fits, pdfs, alpha_eff_pdfs, fits_alpha_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """Plot the central value and the uncertainty of a list of effective
    exponents as a function of x for a given value of Q. If normalize_to
    is given, plot the ratios to the corresponding alpha effective.
    Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the alpha effective or its
    corresponding index in the list, starting from one, or None to plot
    absolute values.
    """
    yield from ExponentBandPlotter(
        fits_alpha_lines, 'alpha', fits, pdfs, alpha_eff_pdfs, 'log', normalize_to, ybottom, ytop)

alpha_eff_fits = collect('alpha_eff', ('fits', 'fitpdfandbasis',))

@figuregen
def plot_alphaEff(
        fits, fits_pdf, alpha_eff_fits, fits_alpha_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """Plot the central value and the uncertainty of a list of effective
    exponents as a function of x for a given value of Q. If normalize_to
    is given, plot the ratios to the corresponding alpha effective.
    Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the alpha effective or its
    corresponding index in the list, starting from one, or None to plot
    absolute values.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.
    """
    return plot_alphaEff_internal(
        fits, fits_pdf, alpha_eff_fits, fits_alpha_lines, normalize_to, ybottom, ytop)

beta_eff_pdfs = collect('beta_eff', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_betaEff_internal(
        fits, pdfs, beta_eff_pdfs, fits_beta_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Same as plot_alphaEff_internal but for beta effective exponent """
    yield from ExponentBandPlotter(
        fits_beta_lines, 'beta', fits, pdfs, beta_eff_pdfs, 'linear', normalize_to, ybottom, ytop)

beta_eff_fits = collect('beta_eff', ('fits', 'fitpdfandbasis',))

@figuregen
def plot_betaEff(
        fits, fits_pdf, beta_eff_fits, fits_beta_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Same as plot_alphaEff but for beta effective exponents """
    return plot_betaEff_internal(
        fits, fits_pdf, beta_eff_fits, fits_beta_lines, normalize_to, ybottom, ytop)

@table
@make_argcheck(check_basis)
def effective_exponents_table_internal(fit: FitSpec, pdf: PDF, *,
                                       x1_alpha: numbers.Real = 1e-6,
                                       x2_alpha: numbers.Real = 1e-3,
                                       x1_beta: numbers.Real = 0.65,
                                       x2_beta: numbers.Real = 0.95,
                                       basis:(str, Basis),
                                       flavours: (list, tuple, type(None)) = None):
    """Returns a table with the effective exponents for the next fit
    By default `x1_alpha = 1e-6`, `x2_alpha = 1e-3`, `x1_beta = 0.65`, and `x2_beta = 0.95`, but
    different values can be specified in the runcard. The values control where the bounds of alpha
    and beta are evaluated:

    alpha_min = singlet/gluon: the 2x68% c.l. lower value evaluated at x=`x1_alpha`
    others  : min(2x68% c.l. lower value evaluated at x=`x1_alpha` and x=`x2_alpha`)
    alpha_max = singlet/gluon: min(2 and the 2x68% c.l. upper value evaluated at x=`x1_alpha`)
                   others    : min(2 and max(2x68% c.l. upper value
                                   evaluated at x=`x1_alpha` and x=`x2_alpha`))
    beta_min  =  max(0 and min(2x68% c.l. lower value evaluated at x=`x1_beta` and x=`x2_beta`))
    beta_max  =  max(2x68% c.l. upper value evaluated at x=`x1_beta` and x=`x2_beta`)
    """

    #Reading from the filter
    with open(fit.path/'filter.yml', 'r') as f:
        filtermap = yaml.safe_load(f)
    previous_exponents = filtermap['fitting']['basis']
    Qmin = pdf.QMin

    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    alpha_effs = alpha_eff(
        pdf, xmin=x1_alpha, xmax=x2_alpha, npoints=2, Q=Qmin, basis=basis, flavours=flavours)
    beta_effs = beta_eff(
        pdf, xmin=x1_beta, xmax=x2_beta, npoints=2, Q=Qmin, basis=basis, flavours=flavours)

    eff_exp_data = []

    alphagrid = alpha_effs.grid_values
    betagrid = beta_effs.grid_values

    alphastats = pdf.stats_class(alphagrid)
    betastats = pdf.stats_class(betagrid)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)

        alpha_cv = np.nanmean(alphagrid, axis=0)
        beta_cv = np.nanmean(betagrid, axis=0)
        #tuple of low and high values repectively
        alpha68 = alphastats.errorbar68()
        beta68 = betastats.errorbar68()

    alpha_sigup = alpha68[1] - alpha_cv
    beta_sigup = beta68[1] - beta_cv
    alpha_sigdown = -alpha68[0] + alpha_cv
    beta_sigdown = -beta68[0] + beta_cv
    flavours_label = []
    runcard_flavours = basis.to_known_elements(
        [ref_fl['fl'] for ref_fl in previous_exponents]).tolist()
    for (j, fl) in enumerate(flavours):

        prev_a_bounds = previous_exponents[runcard_flavours.index(fl)]['smallx']
        prev_b_bounds = previous_exponents[runcard_flavours.index(fl)]['largex']

        # the gluon/singlet case
        if fl in (r'\Sigma', "g"):
            new_alpha_bounds = [
                alpha_cv[j, 0] - 2*alpha_sigdown[j, 0],
                min(2, alpha_cv[j, 0] + 2*alpha_sigup[j, 0])]
        else:
            new_alpha_bounds = [
                min(alpha_cv[j, :] - 2*alpha_sigdown[j, :]),
                min(2, max(alpha_cv[j, :] + 2*alpha_sigup[j, :]))]
        alpha_line = [*prev_a_bounds, *new_alpha_bounds]

        new_beta_bounds = [
            max(0, min(beta_cv[j, :] - 2*beta_sigdown[j, :])),
            max(beta_cv[j, :] + 2*beta_sigup[j, :])]
        beta_line = [*prev_b_bounds, *new_beta_bounds]

        eff_exp_data.extend((alpha_line, beta_line))
        flavours_label.append(f'${fl}$')
    ind = pd.MultiIndex.from_product([flavours_label, [r"$\alpha$", r"$\beta$"]])
    eff_exp_columns = ["prev Min", "prev Max", "next Min", "next Max"]
    df = pd.DataFrame(eff_exp_data, index=ind,
                      columns=eff_exp_columns)
    df.name = pdf.name
    return df


effective_exponents_table = collect(
    'effective_exponents_table_internal', ('fitpdfandbasis',))
fmt = lambda a: float(significant_digits(a, 4))

def next_effective_exponents_yaml(fit: FitSpec, effective_exponents_table):
    """-Returns a table in yaml format called NextEffExps.yaml
       -Prints the yaml table in the report
    using `effective_exponents_table` this provider outputs the yaml runcard to run
    a fit with identical input as the specified `fit` with the t0 and preprocessing iterated.

    This action must be used in a report and should be wrapped in a code block to be formatted
    correctly, for example:

    ```yaml
    {@next_effective_exponents_runcard@}
    ```

    """

    df_effexps = effective_exponents_table[0]
    #Reading from the filter
    with open(fit.path/'filter.yml', 'r') as f:
        filtermap = yaml.load(f, yaml.RoundTripLoader)
    previous_exponents = filtermap['fitting']['basis']
    basis = filtermap['fitting']['fitbasis']
    checked = check_basis(basis, None)
    basis = checked['basis']
    flavours = checked['flavours']

    runcard_flavours = basis.to_known_elements(
        [ref_fl['fl'] for ref_fl in previous_exponents]).tolist()
    for fl in flavours:
        alphas = df_effexps.loc[(f'${fl}$', r'$\alpha$'), ['next Min', 'next Max']].values
        betas = df_effexps.loc[(f'${fl}$', r'$\beta$'), ['next Min', 'next Max']].values
        previous_exponents[runcard_flavours.index(fl)]['smallx'] = [fmt(alpha) for alpha in alphas]
        previous_exponents[runcard_flavours.index(fl)]['largex'] = [fmt(beta) for beta in betas]
    #iterate t0
    filtermap['datacuts']['t0pdfset'] = fit.name
    return yaml.dump(filtermap, Dumper=yaml.RoundTripDumper)
