# -*- coding: utf-8 -*-
"""
Figures for visualizing results
"""
from __future__ import generator_stop

import logging

import numpy as np

from reportengine.figure import figuregen

from validphys.checks import check_scale
from validphys.plots import check_pdf_normalize_to
from validphys.plots import BandPDFPlotter
from validphys.plots import PDFPlotter

from validphys.pdfgrids import xgrid
from validphys.pdfgrids import xplotting_grid

log = logging.getLogger(__name__)

class PreprocessingPlotter(PDFPlotter):
    def __init__(self, exponent, *args,  **kwargs):
        self.exponent = exponent
        super().__init__(*args, **kwargs)

    def get_title(self, parton_name):
        return fr"$\{self.exponent}_e$ for ${parton_name}$ at {self.Q:.1} Gev"

    def get_ylabel(self, parton_name):
        return fr"$\{self.exponent}_e$ for ${parton_name}$"

class ExponentBandPlotter(BandPDFPlotter, PreprocessingPlotter): pass

def alpha_eff(xplotting_grids):
    alpha_grids=[]
    for xplotting_grid in xplotting_grids:
        newvalues = xplotting_grid.grid_values # grid_values[replica][flavour][x]
        xgrid = xplotting_grid.xgrid

        newvalues = -np.log(abs(newvalues/xgrid))/np.log(xgrid)
        #newgrid is like the old grid but with updated values
        newgrid = xplotting_grid._replace(grid_values=newvalues)
        alpha_grids.append(newgrid)

    # xplotting_grid._replace(grid_values=newvalues)
    return alpha_grids #.grid_values

def beta_eff(xplotting_grids):
    beta_grids=[]
    for xplotting_grid in xplotting_grids:
        newvalues = xplotting_grid.grid_values # grid_values[replica][flavour][x]
        xgrid = xplotting_grid.xgrid

        newvalues = np.log(abs(newvalues/xgrid))/np.log(1-xgrid)
        #newgrid is like the old grid but with updated values
        newgrid = xplotting_grid._replace(grid_values=newvalues)
        beta_grids.append(newgrid)

    # xplotting_grid._replace(grid_values=newvalues)
    return beta_grids #.grid_values


@figuregen
@check_pdf_normalize_to
def plot_alphaEff(pdfs, normalize_to:(int,str,type(None))=None,xmin=1e-5,Q=1.65,basis='evolution',flavours=None):

    xpdfs = []
    xmax=0.1
    temp_xgrid = xgrid(xmin, xmax,'log', 200)

    for pdf in pdfs:
        xpdf = xplotting_grid(pdf, Q, xgrid=temp_xgrid, basis=basis,flavours=flavours)
        xpdfs.append(xpdf)

    yield from ExponentBandPlotter('alpha', pdfs, alpha_eff(xpdfs), 'log', normalize_to)

@figuregen
@check_pdf_normalize_to
def plot_betaEff(pdfs, normalize_to:(int,str,type(None))=None,xmin=1e-2,Q=1.65,basis='evolution',flavours=None):

    xpdfs = []
    xmax=0.9
    temp_xgrid = xgrid(xmin, xmax,'linear', 200)

    for pdf in pdfs:
        xpdf = xplotting_grid(pdf, Q, xgrid=temp_xgrid, basis=basis,flavours=flavours)
        xpdfs.append(xpdf)

    yield from ExponentBandPlotter('beta', pdfs, beta_eff(xpdfs), 'linear', normalize_to)
