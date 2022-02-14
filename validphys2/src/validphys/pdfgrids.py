"""
High level providers for PDF and luminosity grids, formatted in such a way
to facilitate plotting and analysis.
"""
from collections import namedtuple
import dataclasses
import numbers
import logging

import numpy as np
import scipy.integrate as integrate

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError, check_positive, check

from validphys.core import PDF, Stats
from validphys.gridvalues import (evaluate_luminosity)
from validphys.pdfbases import (Basis, check_basis)
from validphys.checks import check_pdf_normalize_to, check_xlimits

log = logging.getLogger(__name__)

@make_argcheck
def _check_scale(scale):
    scales = ('linear', 'log')
    if scale not in scales:
        raise CheckError(f'Unrecognized scale {scale}.', scale, scales)

@_check_scale
@check_xlimits
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


@dataclasses.dataclass
class XPlottingGrid:
    """DataClass holding the value of the PDF at the specified
    values of x, Q and flavour.
    The `grid_values` attribute corresponds to a `Stats` instance
    in order to compute statistical estimators in a sensible manner.
    """
    Q: float
    basis: (str, Basis)
    flavours: (list, tuple, type(None))
    xgrid: np.ndarray
    grid_values: Stats
    scale: str

    def copy_grid(self, grid_values=None):
        """Create a copy of the grid with potentially a different set of values"""
        if not isinstance(grid_values, Stats):
            log.warning("XPlottingGrid being called with a numpy grid, should be using Stats instead!")
            grid_values = self.grid_values.__class__(grid_values)
        return dataclasses.replace(self, grid_values=grid_values)


@make_argcheck(check_basis)
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
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']
    if xgrid is None:
        #Call the function that is shadowed
        xgrid = globals()['xgrid']()
    if isinstance(xgrid, tuple) and len(xgrid)==2:
        scale, xgrid = xgrid
    elif isinstance(xgrid, np.ndarray):
        scale = 'unknown'
    else:
        raise TypeError(f"Invalid xgrid {xgrid!r}")
    gv = basis.grid_values(pdf, flavours, xgrid, Q)
    #Eliminante Q axis
    stats_gv = pdf.stats_class(gv.reshape(gv.shape[:-1]))

    res = XPlottingGrid(Q, basis, flavours, xgrid, stats_gv, scale)
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
    stats class, to compute statistics regardless of the error_type.
    """
    s = sqrts*sqrts
    mxs = np.logspace(1, np.log10(sqrts), nbins_m)


    ys  = np.linspace(0 , y_lim, nbins_y)

    y_kinlims = -np.log(mxs/sqrts)
    ys_max = np.searchsorted(ys, y_kinlims)

    # TODO: Write this in something fast
    lpdf = pdf.load()
    nmembers = pdf.get_members()

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

def _default_mxmax(sqrts):
    return sqrts / 3

@make_argcheck
def _check_mx(mxmin, mxmax, sqrts):
    if mxmax is None:
        mxmax = _default_mxmax(sqrts)

    check(
        0 <= mxmin < mxmax <= sqrts,
        ("mxmin and mxmax not consistent: Should be 0 <= mxmin < mxmax <= sqrts, "
        f"but mxmin={mxmin} GeV, mxmax={mxmax} GeV and sqrts={sqrts} GeV."
        ),
    )

@_check_mx
@check_positive("sqrts")
@_check_scale
@check_positive("nbins_m")
def lumigrid1d(
    pdf: PDF,
    lumi_channel,
    sqrts: numbers.Real,
    y_cut: (type(None), numbers.Real) = None,
    nbins_m: int = 50,
    mxmin: numbers.Real = 10,
    mxmax: (type(None), numbers.Real) = None,
    scale="log",
):
    """
    Return the integrated luminosity in a grid of nbins_m points, for the
    values of invariant mass given (proton-proton) collider energy ``sqrts``
    (given in GeV). A rapidity cut on the integration range (if specified)
    is taken into account.

    By default, the grid is sampled logarithmically in mass. The limits are
    given by ``mxmin`` and ``mxmax``, given in GeV. By default ``mxmin`` is 10
    GeV and ``mxmax`` is set based on ``sqrts``.

    The results are computed for all relevant PDF members and wrapped in a
    stats class, to compute statistics regardless of the error_type.
    """
    s = sqrts * sqrts
    if mxmax is None:
        mxmax = _default_mxmax(sqrts)
    if scale == "log":
        mxs = np.logspace(np.log10(mxmin), np.log10(mxmax), nbins_m)
    elif scale == "linear":
        mxs = np.linspace(mxmin, mxmax, nbins_m)
    else:
        raise ValueError("Unknown scale")
    sqrt_taus = (mxs / sqrts)

    # TODO: Write this in something fast
    lpdf = pdf.load()
    nmembers = pdf.get_members()

    weights = np.full(shape=(nmembers, nbins_m), fill_value=np.NaN)

    for im, (mx, sqrt_tau) in enumerate(zip(mxs, sqrt_taus)):
        y_min = -np.log(1/sqrt_tau)
        y_max =  np.log(1/sqrt_tau)

        if y_cut is not None:
            if -y_cut > y_min and  y_cut < y_max:
                y_min = -y_cut
                y_max =  y_cut

        for irep in range(nmembers):
            # Eq.(3) in arXiv:1607.01831
            f = lambda y: evaluate_luminosity(
                lpdf, irep, s, mx,
                sqrt_tau * np.exp(y), sqrt_tau * np.exp(-y),
                lumi_channel
            )
            res = integrate.quad(f, y_min, y_max, epsrel=5e-4, limit=50)[0]

            weights[irep, im] = res

    return Lumi1dGrid(mxs, pdf.stats_class(weights))


lumigrids1d = collect('lumigrid1d', ['lumi_channels'])
pdfs_lumis = collect('lumigrid1d', ('pdfs',))


@check_pdf_normalize_to
def distance_grids(pdfs, xplotting_grids, normalize_to:(int,str,type(None))=None):
    """Return an object containing the value of the distance PDF at the specified values
    of x and flavour.

    The parameter ``normalize_to`` identifies the reference PDF set with respect to the
    distance is computed.

    This method returns distance grids where the relative distance between both PDF
    set is computed. At least one grid will be identical to zero.
    """

    gr2_stats = xplotting_grids[normalize_to].grid_values
    cv2 = gr2_stats.central_value()
    sg2 = gr2_stats.std_error()
    N2 = pdfs[normalize_to].get_members()

    newgrids = []
    for grid, pdf in zip(xplotting_grids, pdfs):

        if pdf == pdfs[normalize_to]:
            # Zero the PDF we are normalizing against
            pdf_zero = pdf.stats_class(np.zeros_like(gr2_stats.data[0]))
            newgrid = grid.copy_grid(grid_values=pdf_zero)
            newgrids.append(newgrid)
            continue

        g_stats = grid.grid_values
        cv1 = g_stats.central_value()
        sg1 = g_stats.std_error()
        N1 = pdf.get_members()

        # the distance definition
        distance = pdf.stats_class(np.sqrt((cv1-cv2)**2/(sg1**2/N1+sg2**2/N2)))

        newgrid = grid.copy_grid(grid_values=distance)
        newgrids.append(newgrid)

    return newgrids


@check_pdf_normalize_to
def variance_distance_grids(pdfs, xplotting_grids, normalize_to:(int,str,type(None))=None):
    """Return an object containing the value of the variance distance PDF at the specified values
    of x and flavour.

    The parameter ``normalize_to`` identifies the reference PDF set with respect to the
    distance is computed.

    This method returns distance grids where the relative distance between both PDF
    set is computed. At least one grid will be identical to zero.
    """

    gr2_stats = xplotting_grids[normalize_to].grid_values
    sg2 = gr2_stats.std_error()
    mo2 = gr2_stats.moment(4)
    N2 = pdfs[normalize_to].get_members()
    s2 = (mo2-(N2-3)/(N2-1)*sg2**4)/N2

    newgrids = []
    for grid, pdf in zip(xplotting_grids, pdfs):

        if pdf == pdfs[normalize_to]:
            # Zero the PDF we are normalizing against
            pdf_zero = pdf.stats_class(np.zeros_like(gr2_stats.data[0]))
            newgrid = grid.copy_grid(grid_values=pdf_zero)
            newgrids.append(newgrid)
            continue

        g_stats = grid.grid_values
        sg1 = g_stats.std_error()
        mo1 = g_stats.moment(4)
        N1 = pdf.get_members()
        s1 = (mo1-(N1-3)/(N1-1)*sg1**4)/N1

        # the distance definition
        variance_distance = pdf.stats_class(np.sqrt((sg1**2-sg2**2)**2/(s1+s2)))

        newgrid = grid.copy_grid(grid_values=variance_distance)
        newgrids.append(newgrid)

    return newgrids
