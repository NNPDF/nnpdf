# -*- coding: utf-8 -*-
"""
Tools for sampling, PDFs.

@author: Zahari Kassabov
"""
from collections import namedtuple

import numpy as np

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError
from NNPDF import _lhapdfset, LHAPDFSet

from validphys.core import PDF

PDG_PARTONS = dict((
                (-5 , r"\bar{b}"),
                (-4 , r"\bar{c}"),
                (-3 , r"\bar{s}"),
                (-2 , r"\bar{u}"),
                (-1 , r"\bar{d}"),
                (0 , r"g"),
                (1 , r"d"),
                (2 , r"u"),
                (3 , r"s"),
                (4 , r"c"),
                (5 , r"b"),
                (22 , r"\gamma"),
                (21 , r"g"),
              )
              )

PDG_ALIASES = {
 '\\bar{b}': -5,
 'bbar'    : -5,
 '\\bar{c}': -4,
 'cbar'    : -4,
 '\\bar{d}': -1,
 'dbar'    : -1,
 '\\bar{s}': -3,
  'sbar'   : -3,
 '\\bar{u}': -2,
  'ubar'   : -2,
 '\\gamma': 22,
 'photon': 22,
 'b': 5,
 'bottom': 5,
 'c': 4,
 'charm': 4,
 'd': 1,
 'down': 1,
 'g': 21,
 'gluon': 21,
 's': 3,
 'strange': 3,
 'u': 2,
 'up': 2,
 }


#Numpy is unhappy about downcasting double to float implicitly, so we have
#to manually make sure all input arrays correspond to NNPDF::real.
FLTYPE = np.int32
REALTYPE = np.double if _lhapdfset.REALDOUBLE else np.float32


DEFAULT_FLARR = (-3,-2,-1,0,1,2,3,4)

def _parse_flarr(flarr):
    out = []
    for elem in flarr:
        msg = "Unknown parton '%s'" % elem
        try:
            num = int(elem)
        except (ValueError, TypeError):
            if elem in PDG_ALIASES:
                out.append(PDG_ALIASES[elem])
            else:
                raise ValueError(msg)
        else:
            if num in PDG_PARTONS:
                out.append(num)
            else:
                raise ValueError(msg)
    return out


def grid_values(pdf:PDF, flmat, xmat, qmat):
    """Returns a 4-dimension array with the PDF values at the input parameters
    for each replica. The return value is indexed as follows:

    grid_values[replica][flavour][x][Q]

    This uses libnnpdf, and therefore follows the convention
    to throw away replica 0 for Monte Carlo ensembles
    (so index 0 corresponds to replica 1). The higher level function
    `central_and_error_grid_values` sets this straight.
    """
    flmat = np.atleast_1d(np.asanyarray(flmat, dtype=FLTYPE))
    xmat, qmat =  (np.atleast_1d(np.asarray(x, dtype=REALTYPE))
                          for x in (xmat,qmat))
    lpdf = pdf.load()
    return lpdf.grid_values(flmat, xmat, qmat)



ScaleSpec = namedtuple('ScaleSpec', ('scale', 'values'))

def _check_xgrid(xgrid):
    #Checking code is always ugly...
    msg = "Incorrect xgrid specification"
    if xgrid in ('log', None):
        xgrid = ('log', np.logspace(-5,0,200, endpoint=False))

    elif xgrid == 'linear':
        xgrid = ('linear', np.linspace(1e-6, 1, 200, endpoint=False))
    elif isinstance(xgrid, list):
        if len(xgrid) in (2,3):
            try:
                cond = xgrid[0] < 0 or xgrid[1] > 1 or xgrid[0] > xgrid[1]
            except Exception as e:
                raise CheckError(msg) from e
            if cond:
                raise CheckError("Incorrect specification: x must be in (0,1) with xmin<=xmax")
        try:
            if len(xgrid) == 2:

                xgrid = ('linear', np.linspace(*xgrid, 200, endpoint=False))
            elif len(xgrid) == 3:
                xgrid = ('linear', np.linspace(*xgrid, endpoint=False))
            else:
                raise CheckError(msg)
        except Exception as e:
            raise CheckError(msg) from e
    else:
        raise CheckError(msg)
    return {'xgrid':xgrid}

def _check_flavours(flavours):
    try:
        return {'flavours': _parse_flarr(flavours)}
    except ValueError as e:
        raise CheckError(e) from e

XPlottingGrid = namedtuple('XPlottingGrid', ('Q', 'flavours', 'xgrid',
                                             'grid_values', 'scale'))


@make_argcheck(_check_xgrid)
@make_argcheck(_check_flavours)
def xplotting_grid(pdf:PDF, Q:(float,int), xgrid=None,
                   flavours:(list, tuple)=DEFAULT_FLARR):
    """Return an object containing the value of the PDF at the specified values
    of x and flavour.

    Q: The PDF scale in GeV.

    The x grid can be specified as follows:
        log: Sample log-spaced points (default).
        linear: Sample linearly-spaced points.
        [xmin,xmax]: linearly sample between xmin and xmax.
        [xmin,xmax,num]: linearly sample num points between xmin and xmax.
    """
    #Make usable outside reportengine
    if isinstance(xgrid, tuple) and len(xgrid)==2:
        scale, xgrid = xgrid
    elif isinstance(xgrid, np.ndarray):
        scale = 'unknown'
    else:
        scale, xgrid = _check_xgrid(xgrid)['xgrid']
    gv = grid_values(pdf, flavours, xgrid, Q)
    #Eliminante Q axis
    #TODO: wrap this in pdf.stats_class?
    gv = gv.reshape(gv.shape[:-1])

    res = XPlottingGrid(Q, flavours, xgrid, gv, scale)
    return res

xplotting_grids = collect(xplotting_grid, ('pdfs',))


LUMI_CHANNELS = {
    'gg': r'gg',
    'gq': r'gq',
    'qqbar': r'q\bar{q}',
    'qq': r'qq',
    'udbar': r'u\bar{d}',
    'dubar': r'd\bar{u}',
}

LUMI_FLAVORS = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]


def _check_channel(channel):
    if channel is None:
        raise CheckError('Channel is None.')

    if channel not in LUMI_CHANNELS:
        raise CheckError('Channel %s not allowed.' % channel)


def evaluate_luminosity(pdf_set: LHAPDFSet, n: int, s: float, mx: float,
                        x1: float, x2: float, channel=None):
    """Returns PDF luminosity at specified values of mx, x1, x2, sqrts**2
    for a given channel.

    pdf_set: The interested PDF set
    s: The center of mass energy GeV.
    mx: The invariant mass bin GeV.
    x1 and x2: The the partonic x1 and x2.
    channel: The channel tag name from LUMI_CHANNELS.
    """

    pdfs = 0
    if channel == 'gg':
        pdfs = pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, 21)
    elif channel == 'gq':
        for i in LUMI_FLAVORS:
            pdfs += pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, 21) \
                    + pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, i)
    elif channel == 'qqbar':
        for i in LUMI_FLAVORS:
            pdfs += pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, -i)
    elif channel == 'qq':
        for i in LUMI_FLAVORS:
            for j in LUMI_FLAVORS:
                pdfs += pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, j)
    elif channel == 'udbar':
        pdfs = pdf_set.xfxQ(x1, mx, n, 2) * pdf_set.xfxQ(x2, mx, n, -1) \
               + pdf_set.xfxQ(x1, mx, n, -1) * pdf_set.xfxQ(x2, mx, n, 2)
    elif channel == 'dubar':
        pdfs = pdf_set.xfxQ(x1, mx, n, 1) * pdf_set.xfxQ(x2, mx, n, -2) \
               + pdf_set.xfxQ(x1, mx, n, -2) * pdf_set.xfxQ(x2, mx, n, 1)

    else:
        raise ValueError("Bad channel")

    return 1.0/s*pdfs/x1/x2


@make_argcheck(_check_channel)
def evaluate_luminosity_grid(pdf: PDF, s: float, mx: float, y: float, channel=None):
    """Returns central value and error of PDF luminosity at the specified values
    of mx, rapidity, sqrts**2 for a given channel.

    pdf: The interested PDF
    s: The center of mass energy GeV.
    mx: The invariant mass bin GeV.
    y: The rapidity value (used to compute the partonic x1 and x2).
    channel: The channel tag name from LUMI_CHANNELS.
    """

    x1 = mx/s**0.5*np.exp(y)
    x2 = mx/s**0.5*np.exp(-y)

    pdf_set = pdf.load()
    members = pdf_set.GetMembers()
    obs = np.zeros(members)

    for n in range(members):
        obs[n] = evaluate_luminosity(pdf_set, n, s, mx, x1, x2, channel)

    uset = pdf.stats_class(obs.reshape(members, 1))
    return uset.central_value(), uset.std_error()
