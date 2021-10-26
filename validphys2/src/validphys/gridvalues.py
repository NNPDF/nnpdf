"""
gridvalues.py

Core functionality needed to obtain a set of values from
LHAPDF. The tools for representing these grids are in pdfgrids.py
(the validphys provider module), and the
basis transformations are in pdfbases.py
"""
import itertools

import numpy as np

from NNPDF import REALDOUBLE, LHAPDFSet

from validphys.core import PDF

#Numpy is unhappy about downcasting double to float implicitly, so we have
#to manually make sure all input arrays correspond to NNPDF::real.
FLTYPE = np.int32
REALTYPE = np.double if REALDOUBLE else np.float32

# Canonical ordering of PDG quark flavour codes
QUARK_FLAVOURS = (-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6)

LUMI_CHANNELS = {
    'gg': r'gg',
    'gq': r'gq',
    'qqbar': r'q\bar{q}',
    'qq': r'qq',
    'ddbar': r'd\bar{d}',
    'uubar': r'u\bar{u}',
    'ssbar': r's\bar{s}',
    'ccbar': r'c\bar{c}',
    'bbbar': r'b\bar{b}',
    'dubar': r'd\bar{u}',
    'udbar': r'u\bar{d}',
    'scbar': r's\bar{c}',
    'csbar': r'c\bar{s}',
    'pp': r'\gamma\gamma',
    'gp': r'g\gamma',
}

QUARK_COMBINATIONS = {
    "ddbar": [1, -1],
    "uubar": [2, -2],
    "ssbar": [3, -3],
    "ccbar": [4, -4],
    "bbbar": [5, -5],
    "dubar": [1, -2],
    "udbar": [2, -1],
    "scbar": [3, -4],
    "csbar": [4, -3],
}

def _grid_values(lpdf, flmat, xmat, qmat):
    """Compute lpdf.grid_values with more forgiving argument types"""
    flmat = np.atleast_1d(np.asanyarray(flmat, dtype=FLTYPE))
    xmat = np.atleast_1d(np.asarray(xmat, dtype=REALTYPE))
    qmat = np.atleast_1d(np.asarray(qmat, dtype=REALTYPE))
    return lpdf.grid_values(flmat, xmat, qmat)

def grid_values(pdf:PDF, flmat, xmat, qmat):
    """
    Evaluate ``x*f(x)`` on a grid of points in flavour, x and Q.

    Parameters
    ----------
    pdf : PDF
        Any PDF set
    flmat : iterable
        A list of PDG IDs corresponding the the LHAPDF flavours in the grid.
    xmat : iterable
        A list of x values
    qmat : iterable
        A list of values in Q, expressed in GeV.

    Returns
    -------
    A 4-dimension array with the PDF values at the input parameters
    for each replica. The return value is indexed as follows::

        grid_values[replica][flavour][x][Q]

    Notes
    ----
    This uses libnnpdf, and therefore follows the convention to throw away
    replica 0 for Monte Carlo ensembles (so index 0 corresponds to replica 1).
    Use ``pdf.grid_values_index`` to index the result properly.

    See Also
    --------
    :py:meth:`validphys.pdfbases.Basis.grid_values` offers a higher level
    interface, allowing to obtain the grid in different bases, and allowing for
    aliases to refer to the flavours.

    Examples
    --------
    Compute the maximum difference across replicas between the u and ubar PDFs
    (times x) for x=0.05 and both Q=10 and Q=100::

        >>> from validphys.loader import Loader
        >>> from validphys.gridvalues import grid_values
        >>> import numpy as np
        >>> gv = grid_values(Loader().check_pdf('NNPDF31_nnlo_as_0118'), [-1, 1], [0.5], [10, 100])
        >>> #Take the difference across the flavour dimension, the max
        >>> #across the replica dimension, and leave the Q dimension untouched.
        >>> np.diff(gv, axis=1).max(axis=0).ravel()
        array([0.07904731, 0.04989902], dtype=float32)
    """
    return _grid_values(pdf.load(), flmat, xmat, qmat)

def central_grid_values(pdf:PDF, flmat, xmat, qmat):
    """Same as :py:func:`grid_values` but it returns only the central values. The
    return value is indexed as::

        grid_values[replica][flavour][x][Q]

    where the first dimension (coresponding to the central member of the PDF set) is
    always one.
    """
    return _grid_values(pdf.load_t0(), flmat, xmat, qmat)


#TODO: Investigate writting these in cython/cffi/numba/...

def evaluate_luminosity(pdf_set: LHAPDFSet, n: int, s:float, mx: float,
                        x1: float, x2: float, channel):
    """Returns PDF luminosity at specified values of mx, x1, x2, sqrts**2
    for a given channel.

    pdf_set: The interested PDF set
    s: The square of the center of mass energy GeV^2.
    mx: The invariant mass bin GeV.
    x1 and x2: The partonic x1 and x2.
    channel: The channel tag name from LUMI_CHANNELS.
    """


    res = 0
    if channel == 'gg':
        res = pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, 21)
    elif channel == 'gq':
        for i in QUARK_FLAVOURS:
            # as in the first of Eq.(4) in arXiv:1607.01831
            res += (pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, 21)
                    + pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, i))
    elif channel == 'gp':
        res = (pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, 22)
               + pdf_set.xfxQ(x1, mx, n, 22) * pdf_set.xfxQ(x2, mx, n, 21))
    elif channel == 'pp':
        res = pdf_set.xfxQ(x1, mx, n, 22) * pdf_set.xfxQ(x2, mx, n, 22)
    elif channel == 'qqbar':
        for i in QUARK_FLAVOURS:
            # as in the third of Eq.(4) in arXiv:1607.01831
            res += pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, -i)
    elif channel == 'qq':
        r1 = []
        r2 = []
        for i in QUARK_FLAVOURS:
            r1.append(pdf_set.xfxQ(x1, mx, n, i))
            r2.append(pdf_set.xfxQ(x2, mx, n, i))

        # as in the second of Eq.(4) in arXiv:1607.01831
        res = sum(a*b for a,b in itertools.product(r1,r2))
    elif channel in QUARK_COMBINATIONS.keys():
        i, j = QUARK_COMBINATIONS[channel]
        res = (pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, j)
               + pdf_set.xfxQ(x1, mx, n, j) * pdf_set.xfxQ(x2, mx, n, i))

    else:
        raise ValueError("Bad channel")

    # The following is equivalent to Eq.(2) in arXiv:1607.01831
    return res/x1/x2/s
