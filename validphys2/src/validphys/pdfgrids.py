"""
High level providers for PDF and luminosity grids, formatted in such a way
to facilitate plotting and analysis.
"""
from collections import namedtuple
import numbers

import numpy as np

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError

from validphys.core import PDF
from validphys.gridvalues import grid_values, evaluate_luminosity
from validphys.pdfbases import PDG_ALIASES, PDG_PARTONS, DEFAULT_FLARR, parse_flarr



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
        return {'flavours': parse_flarr(flavours)}
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
    stats class, to compute statists regardless of the ErrorType.
    """

    s = sqrts**2
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
                res= evaluate_luminosity(lpdf, irep, s,
                    mx, x1, x2, lumi_channel)
                weights[irep, im, iy]  = res


    return Lumi2dGrid(ys, mxs, pdf.stats_class(weights))



lumigrids2d = collect('lumigrid2d', ['lumi_channels'])
