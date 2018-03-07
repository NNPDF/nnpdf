# -*- coding: utf-8 -*-
"""
Tools for computing and plotting effective exponents.
"""
from __future__ import generator_stop

import logging
import warnings
import numpy as np

from reportengine.figure import figuregen
from reportengine.floatformatting import format_number

from validphys.checks import check_scale, CheckError, make_argcheck, check_positive
from validphys.plots import check_pdf_normalize_to
from validphys.plots import BandPDFPlotter
from validphys.plots import PDFPlotter
from validphys.pdfbases import (Basis, check_basis)

import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)

@check_positive('Q')
@pdfgrids._check_limits
@make_argcheck(check_basis)
def alpha_eff(pdfs,xmin=1e-5,xmax=0.1,Q=1.65,basis='evolution',flavours=None):
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
    xGrid = pdfgrids.xgrid(xmin, xmax,'log', 200)

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
def beta_eff(pdfs,xmin=1e-2,xmax=0.9,Q=1.65,basis='evolution',flavours=None):
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
    xGrid = pdfgrids.xgrid(xmin, xmax,'linear', 200)

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
def plot_alphaEff(pdfs, alpha_eff, normalize_to:(int,str,type(None))=None):
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
    yield from ExponentBandPlotter('alpha', pdfs, alpha_eff, 'log', normalize_to)

@figuregen
@check_pdf_normalize_to
def plot_betaEff(pdfs, beta_eff, normalize_to:(int,str,type(None))=None):
    """ Same as plot_alphaEff but for beta effective exponent """
    yield from ExponentBandPlotter('beta', pdfs, beta_eff, 'linear', normalize_to)
