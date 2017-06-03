"""
High level providers for PDF and luminosity grids, formatted in such a way
to facilitate plotting and analysis.
"""
from collections import namedtuple
import numbers

import numpy as np

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError, check_positive

from validphys.core import PDF
from validphys.gridvalues import evaluate_luminosity
from validphys.pdfbases import (parse_flarr, list_bases, Basis, UnknownElement)
import scipy.integrate as integrate


ScaleSpec = namedtuple('ScaleSpec', ('scale', 'values'))


@make_argcheck
def _check_scale(scale):
    scales = ('linear', 'log')
    if scale not in scales:
        raise CheckError(f'Unrecognized scale {scale}.', scale, scales)

@make_argcheck
def _check_limits(xmax, xmin):
    if not (0 <= xmin < xmax <= 1):
        raise CheckError(f'xmin ({xmin}) and xmax ({xmax}) must satisfy \n'
                         '0 <= xmin < xmax <= 1')


@_check_scale
@_check_limits
@check_positive('npoints')
def xgrid(xmin:numbers.Real=1e-5, xmax:numbers.Real=1,
          scale:str='log', npoints:int=200):
    """Return a tuple ``(scale, array)`` where ``scale`` is the input scale
    ("linear" or "log") and ``array`` is generated from the input parameters
    and distributed according to scale."""
    if scale == 'log':
        arr = np.logspace(np.log10(xmin), np.log10(xmax), npoints, endpoint=False)
    elif scale == 'linear':
        arr = np.linspace(xmin, xmax, npoints, endpoint=False)
    return (scale, arr)

def _check_flavours(flavours):
    try:
        return {'flavours': parse_flarr(flavours)}
    except ValueError as e:
        raise CheckError(e) from e

def _check_basis(basis, flavours):
    if isinstance(basis, str):
        bases = list_bases()
        try:
            basis = bases[basis]
        except KeyError:
            raise CheckError(f"Unknown basis '{basis}'", basis, bases)

    if flavours is None:
        flavours = basis.default_elements

    try:
        flavours = basis.to_known_elements(flavours)
    except UnknownElement as e:
        bad = e.args[0]
        raise CheckError(f"Unknown basis element '{bad}'", str(bad),
            alternatives=basis.indexes, display_alternatives='all') from e

    return {'basis':basis, 'flavours':flavours}

XPlottingGrid = namedtuple('XPlottingGrid', ('Q', 'basis', 'flavours', 'xgrid',
                                             'grid_values', 'scale'))


@make_argcheck(_check_basis)
def xplotting_grid(pdf:PDF, Q:(float,int), xgrid=None, basis:(str, Basis)='flavour',
                   flavours:(list, tuple, type(None))=None):
    """Return an object containing the value of the PDF at the specified values
    of x and flavour.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    #Make usable outside reportengine
    checked = _check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']
    if isinstance(xgrid, tuple) and len(xgrid)==2:
        scale, xgrid = xgrid
    elif isinstance(xgrid, np.ndarray):
        scale = 'unknown'
    else:
        raise TypeError(f"Invalid xgrid {xgrid!r}")
    gv = basis.grid_values(pdf, flavours, xgrid, Q)
    #Eliminante Q axis
    #TODO: wrap this in pdf.stats_class?
    gv = gv.reshape(gv.shape[:-1])

    res = XPlottingGrid(Q, basis, flavours, xgrid, gv, scale)
    return res

xplotting_grids = collect(xplotting_grid, ('pdfs',))




Lumi2dGrid = namedtuple('Lumi2dGrid', ['y','m','grid_values'])



def lumigrid2d(pdf:PDF, lumi_channel, sqrts:numbers.Real,
        y_lim:numbers.Real=5, nbins_m:int=100,
        nbins_y:int=50):
    """
    Return the differential luminosity in a grid of (nbins_m x nbins_y)
    points, for the allowed values of invariant mass and rpidity for  given
    (proton-proton) collider energy ``sqrts`` (given in GeV).
    ``y_lim`` specifies the maximum rapidy.

    The grid is sampled linearly in rapidity and logarithmically in mass.

    The results are computed for all relevant PDF members and wrapped in a
    stats class, to compute statistics regardless of the ErrorType.
    """
    s = sqrts*sqrts
    mxs = np.logspace(1, np.log10(sqrts), nbins_m)


    ys  = np.linspace(0 , y_lim, nbins_y)

    y_kinlims = -np.log(mxs/sqrts)
    ys_max = np.searchsorted(ys, y_kinlims)

    # TODO: Write this in something fast
    lpdf = pdf.load()
    nmembers = lpdf.GetMembers()

    weights = np.full(shape=(nmembers, nbins_m, nbins_y), fill_value=np.NaN)

    for irep in range(nmembers):
        for im,mx in enumerate(mxs):
            masked_ys = ys[:ys_max[im]]
            for iy,y in enumerate(masked_ys):
                #TODO: Fill this from lpdf.grid_values?
                x1 = mx/sqrts*np.exp(y)
                x2 = mx/sqrts*np.exp(-y)
                res= evaluate_luminosity(lpdf, irep,
                    s, mx, x1, x2, lumi_channel)
                weights[irep, im, iy]  = res


    return Lumi2dGrid(ys, mxs, pdf.stats_class(weights))


lumigrids2d = collect('lumigrid2d', ['lumi_channels'])


Lumi1dGrid = namedtuple('Lumi1dGrid', ['m','grid_values'])


@check_positive('sqrts')
def lumigrid1d(pdf:PDF, lumi_channel, sqrts:numbers.Real, nbins_m:int=30):
    """
    Return the integrated luminosity in a grid of nbins_m points,
    for the values of invariant mass given (proton-proton) collider
    energy ``sqrts`` (given in GeV).

    The grid is sampled logarithmically in mass.

    The results are computed for all relevant PDF members and wrapped in a
    stats class, to compute statistics regardless of the ErrorType.
    """
    s = sqrts*sqrts
    mxs = np.logspace(1, np.log10(sqrts/10), nbins_m)
    taus = (mxs / sqrts) ** 2

    # TODO: Write this in something fast
    lpdf = pdf.load()
    nmembers = lpdf.GetMembers()

    weights = np.full(shape=(nmembers, nbins_m), fill_value=np.NaN)

    for irep in range(nmembers):
        for im, mx in enumerate(mxs):
            f = lambda x1: evaluate_luminosity(lpdf, irep,
                                               s, mx,
                                               x1, taus[im] / x1,
                                               lumi_channel) / x1
            weights[irep, im] = integrate.quad(f, taus[im], 1.0)[0]

    return Lumi1dGrid(mxs, pdf.stats_class(weights))


lumigrids1d = collect('lumigrid1d', ['lumi_channels'])

