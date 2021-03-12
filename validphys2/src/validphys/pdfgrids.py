"""
High level providers for PDF and luminosity grids, formatted in such a way
to facilitate plotting and analysis.
"""
from collections import namedtuple
import numbers

import numpy as np
import scipy.integrate as integrate

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError, check_positive, check

from validphys.core import PDF
from validphys.gridvalues import (evaluate_luminosity)
from validphys.pdfbases import (Basis, check_basis)
from validphys.checks import check_pdf_normalize_to, check_xlimits

ScaleSpec = namedtuple('ScaleSpec', ('scale', 'values'))


# TODO: store this as a csv?
EXPORT_XGRID = xgrid = np.array([1.00000000000000e-09, 1.29708482343957e-09, 1.68242903474257e-09, 2.18225315420583e-09, 2.83056741739819e-09, 3.67148597892941e-09, 4.76222862935315e-09, 6.17701427376180e-09, 8.01211109898438e-09, 1.03923870607245e-08, 1.34798064073805e-08, 1.74844503691778e-08, 2.26788118881103e-08, 2.94163370300835e-08, 3.81554746595878e-08, 4.94908707232129e-08, 6.41938295708371e-08, 8.32647951986859e-08, 1.08001422993829e-07, 1.40086873081130e-07, 1.81704331793772e-07, 2.35685551545377e-07, 3.05703512595323e-07, 3.96522309841747e-07, 5.14321257236570e-07, 6.67115245136676e-07, 8.65299922973143e-07, 1.12235875241487e-06, 1.45577995547683e-06, 1.88824560514613e-06, 2.44917352454946e-06, 3.17671650028717e-06, 4.12035415232797e-06, 5.34425265752090e-06, 6.93161897806315e-06, 8.99034258238145e-06, 1.16603030112258e-05, 1.51228312288769e-05, 1.96129529349212e-05, 2.54352207134502e-05, 3.29841683435992e-05, 4.27707053972016e-05, 5.54561248105849e-05, 7.18958313632514e-05, 9.31954227979614e-05, 1.20782367731330e-04, 1.56497209466554e-04, 2.02708936328495e-04, 2.62459799331951e-04, 3.39645244168985e-04, 4.39234443000422e-04, 5.67535660104533e-04, 7.32507615725537e-04, 9.44112105452451e-04, 1.21469317686978e-03, 1.55935306118224e-03, 1.99627451141338e-03, 2.54691493736552e-03, 3.23597510213126e-03, 4.09103436509565e-03, 5.14175977083962e-03, 6.41865096062317e-03, 7.95137940306351e-03, 9.76689999624100e-03, 1.18876139251364e-02, 1.43298947643919e-02, 1.71032279460271e-02, 2.02100733925079e-02, 2.36463971369542e-02, 2.74026915728357e-02, 3.14652506132444e-02, 3.58174829282429e-02, 4.04411060163317e-02, 4.53171343973807e-02, 5.04266347950069e-02, 5.57512610084339e-02, 6.12736019390519e-02, 6.69773829498255e-02, 7.28475589986517e-02, 7.88703322292727e-02, 8.50331197801452e-02, 9.13244910278679e-02, 9.77340879783772e-02, 1.04252538208639e-01, 1.10871366547237e-01, 1.17582909372878e-01, 1.24380233801599e-01, 1.31257062945031e-01, 1.38207707707289e-01, 1.45227005135651e-01, 1.52310263065985e-01, 1.59453210652156e-01, 1.66651954293987e-01, 1.73902938455578e-01, 1.81202910873333e-01, 1.88548891679097e-01, 1.95938145999193e-01, 2.03368159629765e-01, 2.10836617429103e-01, 2.18341384106561e-01, 2.25880487124065e-01, 2.33452101459503e-01, 2.41054536011681e-01, 2.48686221452762e-01, 2.56345699358723e-01, 2.64031612468684e-01, 2.71742695942783e-01, 2.79477769504149e-01, 2.87235730364833e-01, 2.95015546847664e-01, 3.02816252626866e-01, 3.10636941519503e-01, 3.18476762768082e-01, 3.26334916761672e-01, 3.34210651149156e-01, 3.42103257303627e-01, 3.50012067101685e-01, 3.57936449985571e-01, 3.65875810279643e-01, 3.73829584735962e-01, 3.81797240286494e-01, 3.89778271981947e-01, 3.97772201099286e-01, 4.05778573402340e-01, 4.13796957540671e-01, 4.21826943574548e-01, 4.29868141614175e-01, 4.37920180563205e-01, 4.45982706956990e-01, 4.54055383887562e-01, 4.62137890007651e-01, 4.70229918607142e-01, 4.78331176755675e-01, 4.86441384506059e-01, 4.94560274153348e-01, 5.02687589545177e-01, 5.10823085439086e-01, 5.18966526903235e-01, 5.27117688756998e-01, 5.35276355048428e-01, 5.43442318565661e-01, 5.51615380379768e-01, 5.59795349416641e-01, 5.67982042055800e-01, 5.76175281754088e-01, 5.84374898692498e-01, 5.92580729444440e-01, 6.00792616663950e-01, 6.09010408792398e-01, 6.17233959782450e-01, 6.25463128838069e-01, 6.33697780169485e-01, 6.41937782762089e-01, 6.50183010158361e-01, 6.58433340251944e-01, 6.66688655093089e-01, 6.74948840704708e-01, 6.83213786908386e-01, 6.91483387159697e-01, 6.99757538392251e-01, 7.08036140869916e-01, 7.16319098046733e-01, 7.24606316434025e-01, 7.32897705474271e-01, 7.41193177421404e-01, 7.49492647227008e-01, 7.57796032432224e-01, 7.66103253064927e-01, 7.74414231541921e-01, 7.82728892575836e-01, 7.91047163086478e-01, 7.99368972116378e-01, 8.07694250750291e-01, 8.16022932038457e-01, 8.24354950923382e-01, 8.32690244169987e-01, 8.41028750298844e-01, 8.49370409522600e-01, 8.57715163684985e-01, 8.66062956202683e-01, 8.74413732009721e-01, 8.82767437504206e-01, 8.91124020497459e-01, 8.99483430165226e-01, 9.07845617001021e-01, 9.16210532771399e-01, 9.24578130473112e-01, 9.32948364292029e-01, 9.41321189563734e-01, 9.49696562735755e-01, 9.58074441331298e-01, 9.66454783914439e-01, 9.74837550056705e-01, 9.83222700304978e-01, 9.91610196150662e-01, 1.00000000000000e+00]).reshape(-1, 1)


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



XPlottingGrid = namedtuple('XPlottingGrid', ('Q', 'basis', 'flavours', 'xgrid',
                                             'grid_values', 'scale'))


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

@make_argcheck
def _check_mx(mxmin, mxmax, sqrts):
    check(
        0 <= mxmin < (mxmax if mxmax is not None else sqrts) <= sqrts,
        "mxmin and mxmax not consistent: Should be 0 <= mxmin < mxmax <= sqrts",
    )

@_check_mx
@check_positive("sqrts")
@check_positive("nbins_m")
def lumigrid1d(
    pdf: PDF,
    lumi_channel,
    sqrts: numbers.Real,
    nbins_m: int = 30,
    mxmin: numbers.Real = 10,
    mxmax: (type(None), numbers.Real) = None,
):
    """
    Return the integrated luminosity in a grid of nbins_m points, for the
    values of invariant mass given (proton-proton) collider energy ``sqrts``
    (given in GeV).

    The grid is sampled logarithmically in mass. The limits are given by
    ``mxmin`` and ``mxmax``, given in GeV. By default ``mxmin`` is 10 GeV and
    ``mxmax`` is set based on ``sqrts``.

    The results are computed for all relevant PDF members and wrapped in a
    stats class, to compute statistics regardless of the ErrorType.
    """
    s = sqrts*sqrts
    if mxmax is None:
        mxmax = sqrts/10
    mxs = np.logspace(np.log10(mxmin), np.log10(mxmax), nbins_m)
    taus = (mxs / sqrts) ** 2

    # TODO: Write this in something fast
    lpdf = pdf.load()
    nmembers = pdf.get_members()

    weights = np.full(shape=(nmembers, nbins_m), fill_value=np.NaN)

    for irep in range(nmembers):
        for im, mx in enumerate(mxs):
            f = lambda x1: evaluate_luminosity(lpdf, irep,
                                               s, mx,
                                               x1, taus[im] / x1,
                                               lumi_channel)
            res = integrate.quad(f, taus[im], 1.0, epsrel=0.05, limit=10)[0]
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

    gr2 = xplotting_grids[normalize_to]
    cv2 = pdfs[normalize_to].stats_class(gr2.grid_values).central_value()
    sg2 = pdfs[normalize_to].stats_class(gr2.grid_values).std_error()
    N2 = gr2.grid_values.shape[0]

    newgrids = list()
    for grid, pdf in zip(xplotting_grids, pdfs):

        if pdf == pdfs[normalize_to]:
            newgrid = grid._replace(grid_values=np.zeros(shape=(grid.grid_values.shape[1], grid.grid_values.shape[2])))
            newgrids.append(newgrid)
            continue

        cv1 = pdf.stats_class(grid.grid_values).central_value()
        sg1 = pdf.stats_class(grid.grid_values).std_error()
        N1 = grid.grid_values.shape[0]

        # the distance definition
        distance = np.sqrt((cv1-cv2)**2/(sg1**2/N1+sg2**2/N2))

        newgrid = grid._replace(grid_values=distance)
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

    gr2 = xplotting_grids[normalize_to]
    sg2 = pdfs[normalize_to].stats_class(gr2.grid_values).std_error()
    mo2 = pdfs[normalize_to].stats_class(gr2.grid_values).moment(4)
    N2 = gr2.grid_values.shape[0]
    s2 = (mo2-(N2-3)/(N2-1)*sg2**4)/N2

    newgrids = list()
    for grid, pdf in zip(xplotting_grids, pdfs):

        if pdf == pdfs[normalize_to]:
            newgrid = grid._replace(grid_values=np.zeros(shape=(grid.grid_values.shape[1], grid.grid_values.shape[2])))
            newgrids.append(newgrid)
            continue

        sg1 = pdf.stats_class(grid.grid_values).std_error()
        mo1 = pdf.stats_class(grid.grid_values).moment(4)
        N1 = grid.grid_values.shape[0]
        s1 = (mo1-(N1-3)/(N1-1)*sg1**4)/N1

        # the distance definition
        variance_distance = np.sqrt((sg1**2-sg2**2)**2/(s1+s2))

        newgrid = grid._replace(grid_values=variance_distance)
        newgrids.append(newgrid)

    return newgrids
