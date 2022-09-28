# -*- coding: utf-8 -*-
"""
Tools for computing and plotting asymptotic exponents.
"""
import logging
import numbers
import warnings

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.figure import figuregen
from reportengine.floatformatting import format_number 
from reportengine.table import table

from validphys.checks import check_positive, check_pdf_normalize_to, make_argcheck, check_xlimits
from validphys.core import PDF
from validphys.pdfbases import check_basis, Basis
from validphys.pdfplots import BandPDFPlotter, PDFPlotter

import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)

@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def alpha_asy(pdf: PDF, *,
              xmin: numbers.Real = 1e-6,
              xmax: numbers.Real = 1e-3,
              npoints: int = 100,
              Q: numbers.Real = 1.65,
              basis: (str, Basis),
              flavours: (list, tuple, type(None)) = None):
    """Returns a list of xplotting_grids containing the value of the asymptotic
    exponent alpha, as defined by the first relationship in Eq. (4) of 
    [arXiv:1604.00024], at the specified value of Q (in GeV), in the interval [xmin, xmax].

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    npoints: the number of sub-intervals in the range [xmin, xmax] on which the 
    derivative is computed.
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
    pdfGrid_values = pdfGrid.grid_values.data
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        dx = np.log(xGrid[1]) - np.log(xGrid[0])
        alphaGrid_values = -np.log(abs(pdfGrid_values))
        alphaGrid_values = np.gradient(alphaGrid_values, dx, axis=2, edge_order=2)
        alphaGrid_values[alphaGrid_values == -np.inf] = np.nan  # when PDF_i =0
        alphaGrid = pdfGrid.copy_grid(grid_values=pdf.stats_class(alphaGrid_values))
        
    return alphaGrid

@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def beta_asy(pdf, *,
             xmin: numbers.Real = 0.6,
             xmax: numbers.Real = 0.9,
             npoints: int = 100,
             Q: numbers.Real = 1.65,
             basis: (str, Basis),
             flavours: (list, tuple, type(None)) = None):
    """Returns a list of xplotting_grids containing the value of the asymptotic
    exponent beta, as defined by the second relationship in Eq. (4) of 
    [arXiv:1604.00024], at the specified value of Q (in GeV), in the interval [xmin, xmax].

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    npoints: the number of sub-intervals in the range [xmin, xmax] on which the 
    derivative is computed.
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
    pdfGrid_values = pdfGrid.grid_values.data
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        dx = xGrid[1] - xGrid[0]
        betaGrid_values = np.log(abs(pdfGrid_values))
        betaGrid_values = (xGrid - 1.) * np.gradient(betaGrid_values, dx, axis=2,edge_order=2)
        betaGrid_values[betaGrid_values == -np.inf] = np.nan  # when PDF_i =0
        betaGrid = pdfGrid.copy_grid(grid_values=pdf.stats_class(betaGrid_values))

    return betaGrid

class AsyExponentBandPlotter(BandPDFPlotter):
    """ Class inheriting from BandPDFPlotter, changing title and ylabel to reflect the asymptotic
    exponent being plotted.
    """

    def __init__(self, exponent, *args,  **kwargs):
        self.exponent = exponent
        super().__init__(*args, **kwargs)

    def get_title(self, parton_name):
        return fr"$\{self.exponent}_a$ for ${parton_name}$ at {format_number(self.Q, 3)} GeV"

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return "Ratio to {}".format(self.normalize_pdf.label)
        else:
            return fr"$\{self.exponent}_a$ for ${parton_name}$"

alpha_asy_pdfs = collect('alpha_asy', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_alpha_asy(
        pdfs, alpha_asy_pdfs, pdfs_alpha_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Plots the alpha asymptotic exponent """
    yield from AsyExponentBandPlotter(
        'alpha',
        pdfs,
        alpha_asy_pdfs,
        'log',
        normalize_to,
        ybottom,
        ytop,)    
        
beta_asy_pdfs = collect('beta_asy', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_beta_asy(
        pdfs, beta_asy_pdfs, pdfs_beta_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Plots the beta asymptotic exponent """
    yield from AsyExponentBandPlotter(
        'beta',
        pdfs,
        beta_asy_pdfs,
        'linear',
        normalize_to,
        ybottom,
        ytop,)

@table
@make_argcheck(check_basis)
def asymptotic_exponents_table(
    pdf: PDF,
    *,
    x_alpha: numbers.Real = 1e-6,
    x_beta: numbers.Real = 0.90,
    Q: numbers.Real = 1.65,   
    basis:(str, Basis),
    flavours: (list, tuple, type(None)) = None,
    npoints=100,):
    """ Returns a table with the values of the asymptotic exponents alpha and beta, as defined 
    in Eq. (4) of [arXiv:1604.00024], at the specified value of x and Q.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    npoints: the number of sub-intervals in the range [xmin, xmax] on which the 
    derivative is computed.
    """

    alpha_a = alpha_asy(
        pdf,
        xmin=x_alpha,
        xmax=1e-3,
        npoints=npoints,
        Q=Q,
        basis=basis,
        flavours=flavours)
    
    beta_a = beta_asy(
        pdf,
        xmin=0.60,
        xmax=x_beta,
        npoints=npoints,
        Q=Q,
        basis=basis,
        flavours=flavours)

    alphastats = alpha_a.grid_values
    betastats = beta_a.grid_values

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)

        alpha_cv = alphastats.central_value()
        beta_cv = betastats.central_value()
        
        alpha_er = alphastats.std_error()
        beta_er = betastats.std_error()

        alpha_68 = alphastats.errorbar68()
        beta_68 = betastats.errorbar68()
    
    flavours_label = []
    asy_exp_mean   = []
    asy_exp_err    = []
    asy_exp_min    = []
    asy_exp_max    = []

    for (j, fl) in enumerate(flavours):
        asy_exp_mean.extend((alpha_cv[j,0],beta_cv[j,-1]))
        asy_exp_err.extend((alpha_er[j,0],beta_er[j,-1]))
        asy_exp_min.extend((alpha_68[0][j,0],beta_68[0][j,-1]))
        asy_exp_max.extend((alpha_68[1][j,0],beta_68[1][j,-1]))
        flavours_label.append(f"${basis.elementlabel(fl)}$")

    asy_exp_data = {"mean": asy_exp_mean, "std": asy_exp_err, "min(68% CL)": asy_exp_min, "max(68% CL)": asy_exp_max}
    ind = pd.MultiIndex.from_product([flavours_label, [r"$\alpha$", r"$\beta$"]])

    df = pd.DataFrame(asy_exp_data, index=ind, columns=["mean","std","min(68% CL)","max(68% CL)"])
    return df
