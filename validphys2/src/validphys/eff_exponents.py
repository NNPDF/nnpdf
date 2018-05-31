# -*- coding: utf-8 -*-
"""
Tools for computing and plotting effective exponents.
"""
from __future__ import generator_stop

import logging
import warnings
import numpy as np
import pandas as pd

from reportengine.figure import figuregen
from reportengine.table  import table
from reportengine.floatformatting import format_number

from validphys.checks import check_scale, CheckError, make_argcheck, check_positive, check_pdf_normalize_to
from validphys.plots import BandPDFPlotter
from validphys.plots import PDFPlotter
from validphys.pdfbases import (Basis, check_basis)

import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)

@check_positive('Q')
@pdfgrids._check_limits
@make_argcheck(check_basis)
def alpha_eff(pdfs,xmin=1e-5,xmax=1e-3,Nstep=200,Q=1.65,basis='evolution',flavours=None):
    """Return a list of xplotting_grids containing the value of the effective
    exponent alpha at the specified values of x and flavour.
    alpha is relevant at small x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    alphaGrids=[]
    xGrid = pdfgrids.xgrid(xmin, xmax,'log', Nstep)

    for pdf in pdfs:
        pdfGrid = pdfgrids.xplotting_grid(pdf, Q, xgrid=xGrid, basis=basis,flavours=flavours)
        pdfGrid_values = pdfGrid.grid_values
        xGrid = pdfGrid.xgrid #NOTE: without this I get "setting an array element with a sequence"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            alphaGrid_values = -np.log(abs(pdfGrid_values/xGrid))/np.log(xGrid)
        alphaGrid = pdfGrid._replace(grid_values=alphaGrid_values)
        alphaGrids.append(alphaGrid)

    return alphaGrids #.grid_values


@check_positive('Q')
@pdfgrids._check_limits
@make_argcheck(check_basis)
def beta_eff(pdfs,xmin=0.5,xmax=0.9,Nstep=200,Q=1.65,basis='evolution',flavours=None):
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

    betaGrids=[]
    xGrid = pdfgrids.xgrid(xmin, xmax,'linear', Nstep)

    for pdf in pdfs:
        pdfGrid = pdfgrids.xplotting_grid(pdf, Q, xgrid=xGrid, basis=basis,flavours=flavours)
        pdfGrid_values = pdfGrid.grid_values
        xGrid = pdfGrid.xgrid #NOTE: without this I get "setting an array element with a sequence"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            betaGrid_values = np.log(abs(pdfGrid_values/xGrid))/np.log(1-xGrid)
        betaGrid = pdfGrid._replace(grid_values=betaGrid_values)
        betaGrids.append(betaGrid)

    return betaGrids #.grid_values

class PreprocessingPlotter(PDFPlotter):
    """ Class inherenting from BandPDFPlotter, has the same functionality
    but with overloaded title and ylabel to take into account the effective
    exponents names.
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

class ExponentBandPlotter(BandPDFPlotter, PreprocessingPlotter): pass

@figuregen
@check_pdf_normalize_to
def plot_alphaEff(pdfs, alpha_eff, normalize_to:(int,str,type(None))=None, ymin=None, ymax=None):
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
    yield from ExponentBandPlotter('alpha', pdfs, alpha_eff, 'log', normalize_to, ymin, ymax)

@figuregen
@check_pdf_normalize_to
def plot_betaEff(pdfs, beta_eff, normalize_to:(int,str,type(None))=None, ymin=None, ymax=None):
    """ Same as plot_alphaEff but for beta effective exponent """
    yield from ExponentBandPlotter('beta', pdfs, beta_eff, 'linear', normalize_to, ymin, ymax)

@table
def effective_exponents_table(pdfs,xmin_alpha=1e-6,xmax_alpha=1e-3,xmin_beta=0.65,xmax_beta=0.95,Q=1.65,basis='evolution',flavours=None):
    """Return a table with the effective exponents for the next fit"""
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    xGridmin_alpha = pdfgrids.xgrid(xmin_alpha, xmin_alpha,'log', 1)
    xGridmax_alpha = pdfgrids.xgrid(xmax_alpha, xmax_alpha,'log', 1)
    xGridmin_beta = pdfgrids.xgrid(xmin_beta, xmin_beta,'log', 1)
    xGridmax_beta = pdfgrids.xgrid(xmax_beta, xmax_beta,'log', 1)

    alphamin_grids=alpha_eff(pdfs,xmin=xmin_alpha,xmax=xmin_alpha,Nstep=1,Q=1.65,basis='evolution',flavours=None)
    alphamax_grids=alpha_eff(pdfs,xmin=xmax_alpha,xmax=xmax_alpha,Nstep=1,Q=1.65,basis='evolution',flavours=None)
    betamin_grids=alpha_eff(pdfs,xmin=xmin_beta,xmax=xmin_beta,Nstep=1,Q=1.65,basis='evolution',flavours=None)
    betamax_grids=alpha_eff(pdfs,xmin=xmax_beta,xmax=xmax_beta,Nstep=1,Q=1.65,basis='evolution',flavours=None)

    eff_exp_data=[]


    for pdf in pdfs:
        i=0
        alphamin_grid_values = alphamin_grids[i].grid_values
        alphamax_grid_values = alphamax_grids[i].grid_values
        betamin_grid_values = betamin_grids[i].grid_values
        betamax_grid_values = betamax_grids[i].grid_values

        alphamin_stats = pdfs[i].stats_class(alphamin_grid_values)
        alphamax_stats = pdfs[i].stats_class(alphamax_grid_values)
        betamin_stats = pdfs[i].stats_class(betamin_grid_values)
        betamax_stats = pdfs[i].stats_class(betamax_grid_values)

        i=i+1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            alphamin_err68down, alphamin_err68up = alphamin_stats.errorbar68()
            alphamax_err68down, alphamax_err68up = alphamax_stats.errorbar68()
            betamin_err68down, betamin_err68up = betamin_stats.errorbar68()
            betamax_err68down, betamax_err68up = betamax_stats.errorbar68()

        #Defining the bounds
        # alpha_min = singlet/gluon: the 2x68% c.l. lower value evaluated at x=1e-6.
        #                  others  : min(2x68% c.l. lower value evaluated at x=1e-6 and x=1e-3)
        # alpha_max = singlet/gluon: min(2 and the 2x68% c.l. upper value evaluated at x=1e-6)
        #                others    : min(2 and max(2x68% c.l. upper value evaluated at x=1e-6 and x=1e-3))
        # beta_min  =  max(0 and min(2x68% c.l. lower value evaluated at x=0.65 and x=0.95))
        # beta_max  =  max(2x68% c.l. upper value evaluated at x=0.65 and x=0.95)
        alpha_line=[]
        beta_line=[]

        if i is 0 or 1: #the gluon/singlet case
            alpha_line.append("alpha")
            alpha_line.append(alphamin_err68down)
            alpha_line.append(min(2,alphamin_err68up))
        else:
            alpha_line.append("alpha")
            alpha_line.append(min(alphamin_err68down,alphamax_err68down))
            alpha_line.append(min(2,max(alphamin_err68up,alphamax_err68up)))

        beta_line.append("beta")
        beta_line.append(max(0,min(betamin_err68down,betamax_err68down)))
        beta_line.append(max(betamin_err68up,betamax_err68up))
        eff_exp_data.append(alpha_line)
        eff_exp_data.append(beta_line)

    flavours_label=[]
    for fl in flavours:
        flavours_label.append(f'${basis.elementlabel(fl)}$')
        flavours_label.append("")

    eff_exp_columns=["Effective exponent","Min","Max"]

    return pd.DataFrame(np.random.randn(16, 3),index=flavours_label,columns=eff_exp_columns)
