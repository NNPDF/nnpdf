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

log = logging.getLogger(__name__)

def alpha_eff(xplotting_grid):
    newvalues = xplotting_grid.grid_values # grid_values[replica][flavour][x]
    xgrid = xplotting_grid.xgrid

    newvalues = -np.log(abs(newvalues/xgrid))/np.log(xgrid)
    #newgrid is like the old grid but with updated values
    newgrid = type(xplotting_grid)(**{**xplotting_grid._asdict(),
                        'grid_values':newvalues})
    return newgrid #.grid_values

def beta_eff(xplotting_grid):
    newvalues = xplotting_grid.grid_values # grid_values[replica][flavour][x]
    xgrid = xplotting_grid.xgrid

    newvalues = np.log(abs(newvalues/xgrid))/np.log(1-xgrid)
    #newgrid is like the old grid but with updated values
    newgrid = type(xplotting_grid)(**{**xplotting_grid._asdict(),
                        'grid_values':newvalues})
    return newgrid #.grid_values


@figuregen
@check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_alphaEff(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):

    yield from BandPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to, alpha_eff, 'alpha effective')

@figuregen
@check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_betaEff(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):

    yield from BandPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to, beta_eff, 'beta effective')
