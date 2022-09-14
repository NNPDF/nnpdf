"""
pdfbases.py

This holds the concrete labels data relative to the PDF bases,
as declaratively as possible.
"""
import inspect
import functools
import abc

import numpy as np

from reportengine.checks import CheckError

from validphys.gridvalues import grid_values, central_grid_values


#This mapping maps the keys passed to LHAPDF (PDG codes) to nice LaTeX labels.
PDG_PARTONS = dict((
        (-6,  r'\bar{t}'),
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
        (6 , r"t"),
        (22 , r"\gamma"),
        (21 , r"g"),
    ))

PIDS_DICT = {
        -6: "TBAR",
        -5: "BBAR",
        -4: "CBAR",
        -3: "SBAR",
        -2: "UBAR",
        -1: "DBAR",
        21: "GLUON",
        1: "D",
        2: "U",
        3: "S",
        4: "C",
        5: "B",
        6: "T",
        22: "PHT",
    }

# Canonical ordering of PDG codes (so flavour basis)
ALL_FLAVOURS = (-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6, 22)
DEFAULT_FLARR = (-3,-2,-1,0,1,2,3,4)

def pdg_id_to_canonical_index(flindex):
    """Given an LHAPDF id, return its index in the ALL_FLAVOURS list."""
    if flindex == 0:
        return ALL_FLAVOURS.index(21)
    return ALL_FLAVOURS.index(flindex)

def list_bases():
    """ List available PDF bases """
    import validphys.pdfbases as thismodule
    return dict(inspect.getmembers(thismodule, lambda x: isinstance(x, thismodule.Basis)))

def check_basis(basis, flavours):
    """
        Check to verify a given basis and set of flavours.
        Returns a dictionary with the relevant instance of the basis
        class and flavour specification
    """

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

#These are various ways to refering to the PDG partons with nicer text labels.
PDG_ALIASES = {
 r'\bar{t}': -6,
 'tbar'    : -6,
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
 't': 6,
 'top': 6,
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
 0 : 21,
 }

def parse_flarr(flarr):
    """Parse a free form list into a list of PDG parton indexes
    (that may contain indexes or values from `PDF_ALIASES`)"""
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

class UnknownElement(KeyError):
    pass


class Basis(abc.ABC):
    """A Basis maps a set of PDF flavours (typically as given by
    :ref:`LHAPDF <lhapdf>`) to functions thereof. This abstract class provides
    functionalities to manage labels (used for plotting) and defaults, while
    the concrete implementation of the transformations is handled by the
    subclasses (by implementing the
    :py:meth:`validphys.pdfbases.Basis.apply_grid_values` method). The high
    level :py:meth:`validphys.pdfbases.Basis.grid_values` and
    :py:meth:`validphys.pdfbases.Basis.central_grid_values` methods then provide
    convenient functionality to work with transformations.

    Attributes
    ----------
    labels: list
        A list of strings representing the labels of each possible
        transformation, in order.
    aliases: dict, optional
        A mapping from strings to labels appearing in ``labels``, specifying
        equivalent ways to enter elements in the user interface.
    default_elements: list, optional
        A list of the labels to be computed by default when no subset of
        elements is specified. If not given it is assumed to be the same as
        ``labels``.
    element_representations: dict, optional
        A mapping from strings to labels indicating the preferred string
        representation of the provided elements (to be used in plotting). If
        this parameter is not given or the element is not in the mapping, the
        label itself is used. It may be convenient to set this when heavy use
        of LaTeX is desired.
    """

    def __init__(self, labels, *, aliases=None,
            default_elements=None, element_representations=None):

        self.labels = labels

        #self._known_flavours = ALL_FLAVOURS[]
        if default_elements is None:
            default_elements = labels
        self.default_elements = default_elements
        if element_representations is None:
            element_representations = {}
        self.element_representations = element_representations

        indexes = {lb:i for i,lb in enumerate(labels)}
        if aliases:
             indexes.update({alias:indexes[alias_label] for alias,alias_label in aliases.items()})
        else:
            aliases = []

        self.aliases = aliases

        self.indexes = indexes

    def elementlabel(self, element):
        """Return the printable representation of a given element of this
        basis."""
        if element in self.aliases:
            element = self.aliases[element]
        if element in self.element_representations:
            return self.element_representations[element]
        elif element in self.labels:
            return element
        raise UnknownElement(element)

    def has_element(self, element):
        """ Return true if basis has knowledge of the given element """
        try:
            self.elementlabel(element)
            return True
        except UnknownElement:
            return False

    def _to_indexes(self, basis_arr):
        """Convert a list of elements of the basis to indexes of the
        (rows of the) transformation matrix."""
        return [self.indexes[k] for k in basis_arr]

    def to_known_elements(self, vmat):
        """Transform the list of aliases into an array of known labels. Raise
        `UnknownElement` on failure."""
        try:
            return np.asanyarray(self.labels)[self._to_indexes(vmat)]
        except KeyError as e:
            raise UnknownElement(*e.args) from e

    @abc.abstractmethod
    def apply_grid_values(self, func, vmat, xmat, qmat):
        """Abstract method to implement basis transformations. It outsources the
        filling of the grid in the flavour basis to ``func`` and implements the
        transformation from the flavour basis to the basis. Methods like
        :py:meth:`validphys.pdfbases.Basis.grid_values` and
        :py:meth:`validphys.pdfbases.Basis.central_grid_values` are derived
        from this method by selecting the appropriate ``func``.

        It should return an array indexed as

            grid_values[N][flavour][x][Q]

        Parameters
        ----------
        func: callable
            A function that fills the grid defined by the rest of the input
            with elements in the flavour basis.
        vmat: iterable
            A list of flavour aliases valid for the basis.
        xmat: iterable
            A list of x values
        qmat: iterable
            A list of values in Q, expressed in GeV.

        """
        ...

    def grid_values(self, pdf, vmat, xmat, qmat):
        """Like :py:func:`validphys.gridvalues.grid_values`, but taking  and
        returning `vmat` in terms of the vectors in this base.

        Parameters
        ----------
        pdf: PDF
            Any PDF set
        vmat: iterable
            A list of flavour aliases valid for the basis.
        xmat: iterable
            A list of x values
        qmat: iterable
            A list of values in Q, expressed in GeV.

        Returns
        -------
        grid: np.ndarray
            A 4-dimension array with the PDF values at the input parameters
            for each replica. The return value is indexed as follows:

                grid_values[replica][flavour][x][Q]

        Examples
        --------
        Compute the median ratio over replicas between singlet and gluon for a
        fixed point in x and a range of values in Q::

            >>> import numpy as np
            >>> from validphys.loader import Loader
            >>> from validphys.pdfbases import evolution
            >>> gv = evolution.grid_values(Loader().check_pdf("NNPDF31_nnlo_as_0118"), ["singlet", "gluon"], [0.01], [2,20,200])
            >>> np.median(gv[:,0,...]/gv[:,1,...], axis=0)
            array([[0.56694959, 0.53782002, 0.60348812]])
        """
        func = functools.partial(grid_values, pdf)
        return self.apply_grid_values(func, vmat, xmat, qmat)

    def central_grid_values(self, pdf, vmat, xmat, qmat):
        """Same as :py:meth:`Basis.grid_values` but returning information on
        the central member of the PDF set."""
        func = functools.partial(central_grid_values, pdf)
        return self.apply_grid_values(func, vmat, xmat, qmat)

class LinearBasis(Basis):
    """A basis that implements a linear transformation of flavours.

    Attributes
    ----------
    from_flavour_mat: np.ndarray
        A matrix that rotates the flavour basis into this basis.
    """

    def __init__(self, labels, from_flavour_mat, *args, **kwargs):
        """
        `flavour_representantion` is a mapping with the printable strings of
        the elements (in case it doesn't match `labels`).
        """
        self.from_flavour_mat = from_flavour_mat
        super().__init__(labels, *args, **kwargs)


    """
    NOTE: At the moment, we don't need the inverse functionality, namely
    transforming from basis to flavour. But if we do, it should be computed
    using arbitrary precision inverses with simpy (simpy.Matrix.pinv).

    Of course, other possible strategies invclude writing the inverse
    transformations by hand or using the SVD but not the inverse directly.

    instead of numpy. With  numpy it's too difficult to get exact cancelations
    and the like. In particular, make sure you pass a test like:

    from hypothesis import given
    from hypothesis.strategies import floats
    import hypothesis.extra.numpy as npmats
    import numpy as np

    from validphys.pdfbases import ALL_FLAVOURS, evolution

    flmats = npmats.arrays(float, len(ALL_FLAVOURS), floats(allow_nan=False,
        allow_infinity=False))

    @given(inp=flamts)
    def test_evolution_transformation(inp):
        assert np.allclose(evolution.to_flavour(evolution.from_flavour(inp)),inp,
                       atol=1e-5, rtol=1e-2)

    for some reasonable tolerances.


    In summary, the following code works but with not so good precision:

    ___________________________

    def to_flavour(self, basis_arr):
        return self.to_flavour_mat @ basis_arr



    @property
    def from_flavour_mat(self):
        if self._from_flavour_mat is None:
            mat = la.pinv(self.to_flavour_mat, rcond=1e-3)
            #TODO: This is a rather ugly hack.  sympy can do pseudoinverses
            #with exact arithmetic so maybe we should just use that if we add
            #it as a dependency. Or perhaps
            #define the inverse transofrms by hand.
            mat[np.abs(mat)<1e-15]=0
            self._from_flavour_mat = mat

        return self._from_flavour_mat

    def from_flavour(self, flavour_arr):
        return self.from_flavour_mat @ flavour_arr

    def to_other_base(self, basis_arr, other_base):
        return other_base.from_flavour_mat @ (self.to_flavour_mat @ basis_arr)

    _____________

    There is also the question of the interface, since we almost never want
    all flavours.

    """




    def _flaray_from_flindexes(self, flinds):
        """Convert a list of flavor basis indexes to PDG codes to pass to
        LHAPDF"""
        return np.asarray(ALL_FLAVOURS)[flinds]

    def _flmask_from_base_indexes(self, inds):
        """Return the flavour indexes of the transformation matrix
        (i.e. columns) that are needed to compute `inds` (indexes of the rows,
        possibly obtained with `_to_indexes`)."""
        return np.count_nonzero(self.from_flavour_mat[inds, :], axis=0) > 0


    def apply_grid_values(self, func, vmat, xmat, qmat):

        #Indexes in the transformation from the "free form" inpt
        inds = self._to_indexes(vmat)

        #Indexes in flavour basis required to compute these elements
        flinds = self._flmask_from_base_indexes(inds)

        #The submatrix of the transformation matrix
        index = np.ix_(inds, flinds)
        transformation = self.from_flavour_mat[index]

        #The PDG codes for LHAPDF
        flmat = self._flaray_from_flindexes(flinds)

        gv = func(flmat, xmat, qmat)
        #Matrix product along the flavour axis
        rotated_gv = np.einsum('bc,acde->abde', transformation, gv)
        return rotated_gv



    @classmethod
    def from_mapping(cls, mapping, *, aliases=None, default_elements=None):
        """Construct a basus from a mapping of the form
        ``{label:{pdf_flavour:coefficient}}``."""
        arr = np.zeros(shape=(len(mapping), len(ALL_FLAVOURS)))
        labels = tuple(mapping)
        for i, coefs in enumerate(mapping.values()):
            indexes = [pdg_id_to_canonical_index(val) for val in parse_flarr(coefs.keys())]
            values = coefs.values()
            arr[i, indexes] = list(values)
        return cls(labels, arr, aliases=aliases, default_elements=default_elements)


class ScalarFunctionTransformation(Basis):
    """A basis that transforms the flavour basis into a single element given by
    ``transform_func``.

    Optional keyword arguments are passed to the constructor of
    :py:class:`validphys.pdfbases.Basis`.

    Attributes
    ----------
    transform_func: callable
        A callable with the signature ``transform_func(func, xmat, qmat)``
        that fills the grid in :math:`x` and :math:`Q`  using ``func``
        and returns a grid with a single basis element.

    """
    def __init__(self, transform_func, *args, **kwargs):
        self.transform_func = transform_func
        super().__init__(*args, **kwargs)

    def apply_grid_values(self, func, vmat, xmat, qmat):
        return self.transform_func(func, xmat, qmat)


def scalar_function_transformation(label, *args, **kwargs):
    """Convenience decorator factory to produce a
    :py:class:`validphys.pdfbases.ScalarFunctionTransformation` basis from a
    function.

    Parameters
    ----------
    label: str
        The single label of the element produced by the function transformation.

    Notes
    -----
    Optional keyword arguments are passed to the constructor of
    :py:class:`validphys.pdfbases.ScalarFunctionTransformation`.

    Returns
    -------
    decorator: callable
       A decorator that can be applied to a suitable transformation function.
    """
    def f_(transform_func):
        return ScalarFunctionTransformation(transform_func, [label], *args, **kwargs)

    return f_


flavour = LinearBasis(ALL_FLAVOURS, np.eye(len(ALL_FLAVOURS)), aliases=PDG_ALIASES,
    default_elements = DEFAULT_FLARR, element_representations=PDG_PARTONS
)

#dicts are oredered in python 3.6+... code shouldn't vreak if they aren't
#though
#see Eqs.(56),(57) https://arxiv.org/pdf/0808.1231.pdf for evolution basis definition
evolution = LinearBasis.from_mapping({
    r'\Sigma'  : {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1, 'c': 1, 'cbar': 1 ,'b':1, 'bbar': 1, 't': 1, 'tbar': 1},
    'V'        : {'u': 1, 'ubar':-1, 'd': 1, 'dbar':-1, 's': 1, 'sbar':-1, 'c': 1, 'cbar':-1 ,'b':1, 'bbar':-1, 't': 1, 'tbar':-1},

    'T3'       : {'u': 1, 'ubar': 1, 'd':-1, 'dbar':-1},
    'V3'       : {'u': 1, 'ubar':-1, 'd':-1, 'dbar': 1},

    'T8'       : {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's':-2, 'sbar':-2},
    'V8'       : {'u': 1, 'ubar':-1, 'd': 1, 'dbar':-1, 's':-2, 'sbar':+2},

    'T15'      : {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1, 'c':-3, 'cbar':-3},
    'V15'      : {'u': 1, 'ubar':-1, 'd': 1, 'dbar':-1, 's': 1, 'sbar':-1, 'c':-3, 'cbar':+3},


    'T24'      : {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1, 'c': 1, 'cbar': 1, 'b':-4, 'bbar':-4},
    'V24'      : {'u': 1, 'ubar':-1, 'd': 1, 'dbar':-1, 's': 1, 'sbar':-1, 'c': 1, 'cbar':-1, 'b':-4, 'bbar':+4},

    'T35'      : {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1, 'c': 1, 'cbar': 1, 'b': 1, 'bbar': 1, 't':-5, 'tbar':-5},
    'V35'      : {'u': 1, 'ubar':-1, 'd': 1, 'dbar':-1, 's': 1, 'sbar':-1, 'c': 1, 'cbar':-1, 'b': 1, 'bbar':-1, 't':-5, 'tbar':+5},

    'g'        : {'g':1},
    'photon'   : {'photon':1},
    },
    aliases = {'gluon':'g', 'singlet': r'\Sigma', 'sng': r'\Sigma', 'sigma': r'\Sigma',
               'v': 'V', 'v3': 'V3', 'v8': 'V8', 't3': 'T3', 't8': 'T8', 't15': 'T15'},
    default_elements=(r'\Sigma', 'V', 'T3', 'V3', 'T8', 'V8', 'T15', 'gluon', )
)

EVOL = evolution

PDF4LHC20 = LinearBasis.from_mapping({
        r'\Sigma': {
            'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1,
            'c': 1, 'cbar': 1, 'b': 1, 'bbar': 1, 't': 1, 'tbar': 1},
        'g': {'g': 1},

        'V': {
            'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': 1, 'sbar': -1,
            'c': 1, 'cbar': -1, 'b': 1, 'bbar': -1, 't': 1, 'tbar': -1},

        'V3': {'u': 1, 'ubar': -1, 'd': -1, 'dbar': 1},

        'T3': {'u': 1, 'ubar': 1, 'd': -1, 'dbar': -1},
        'T8': {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': -2, 'sbar': -2},
    
        'photon': {'photon': 1},
    },
    aliases = {'gluon':'g', 'singlet': r'\Sigma', 'sng': r'\Sigma', 'sigma': r'\Sigma',
               'v': 'V', 'v3': 'V3', 't3': 'T3', 't8': 'T8'},
    default_elements=(r'\Sigma', 'gluon', 'V', 'V3', 'T3', 'T8',  ))

NN31IC = LinearBasis.from_mapping(
    {
        r'\Sigma': {
            'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1,
            'c': 1, 'cbar': 1, 'b': 1, 'bbar': 1, 't': 1, 'tbar': 1},
        'g': {'g': 1},

        'V': {
            'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': 1, 'sbar': -1,
            'c': 1, 'cbar': -1, 'b': 1, 'bbar': -1, 't': 1, 'tbar': -1},

        'V3': {'u': 1, 'ubar': -1, 'd': -1, 'dbar': 1},
        'V8': {'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': -2, 'sbar': +2},

        'T3': {'u': 1, 'ubar': 1, 'd': -1, 'dbar': -1},
        'T8': {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': -2, 'sbar': -2},

        r'c^+': {'c': 1, 'cbar': 1},

        'photon': {'photon': 1},
    },
    aliases={
        'gluon': 'g', 'singlet': r'\Sigma', 'sng': r'\Sigma', 'sigma': r'\Sigma', 'cp': r'c^+',
        'v': 'V', 'v3': 'V3', 'v8': 'V8', 't3': 'T3', 't8': 'T8'},
    default_elements=(r'\Sigma', 'gluon', 'V', 'V3', 'V8', 'T3', 'T8', r'c^+', ))

NN31PC = LinearBasis.from_mapping(
    {
        r'\Sigma': {
            'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1,
            'c': 1, 'cbar': 1, 'b': 1, 'bbar': 1, 't': 1, 'tbar': 1},
        'g': {'g': 1},

        'V': {
            'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': 1, 'sbar': -1,
            'c': 1, 'cbar': -1, 'b': 1, 'bbar': -1, 't': 1, 'tbar': -1},

        'V3': {'u': 1, 'ubar': -1, 'd': -1, 'dbar': 1},
        'V8': {'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': -2, 'sbar': +2},

        'T3': {'u': 1, 'ubar': 1, 'd': -1, 'dbar': -1},
        'T8': {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': -2, 'sbar': -2},
        'photon': {'photon': 1},
    },
    aliases={
        'gluon': 'g', 'singlet': r'\Sigma', 'sng': r'\Sigma', 'sigma': r'\Sigma',
        'v': 'V', 'v3': 'V3', 'v8': 'V8', 't3': 'T3', 't8': 'T8'},
    default_elements=(r'\Sigma', 'gluon', 'V', 'V3', 'V8', 'T3', 'T8'))

FLAVOUR = LinearBasis.from_mapping(
    {
        'u': {'u': 1},
        'ubar': {'ubar': 1},
        'd': {'d': 1},
        'dbar': {'dbar': 1},
        's': {'s': 1},
        'sbar': {'sbar': 1},
        'c': {'c': 1},
        'g': {'g': 1},
    },
    default_elements=('u', 'ubar', 'd', 'dbar', 's', 'sbar', 'c', 'g', ))

pdg = LinearBasis.from_mapping({
'g/10': {'g':0.1},
'u_{v}': {'u':1, 'ubar':-1},
'd_{v}': {'d':1, 'dbar': -1},
's': {'s':1},
r'\bar{u}': {'ubar':1},
r'\bar{d}': {'dbar':1},
'c': {'c':1},
})


@scalar_function_transformation(label="u/d")
def ud_ratio(func, xmat, qmat):
    gv = func([2, 1], xmat, qmat)
    num = gv[:, [0], ...]
    den = gv[:, [1], ...]
    return num / den

@scalar_function_transformation(label="d/u")
def du_ratio(func, xmat, qmat):
    gv = func([1, 2], xmat, qmat)
    num = gv[:, [0], ...]
    den = gv[:, [1], ...]
    return num / den

@scalar_function_transformation(label=r"\bar{d}/\bar{u}")
def dbarubar_ratio(func, xmat, qmat):
    gv = func([-1, -2], xmat, qmat)
    num = gv[:, [0], ...]
    den = gv[:, [1], ...]
    return num / den

  
@scalar_function_transformation(label="Rs", element_representations={"Rs": "R_{s}"})
def strange_fraction(func, xmat, qmat):
    gv = func([-3, 3, -2, -1], xmat, qmat)
    sbar, s, ubar, dbar = (gv[:, [i], ...] for i in range(4))
    return (sbar + s) / (ubar + dbar)

  
def fitbasis_to_NN31IC(flav_info, fitbasis):
    """Return a rotation matrix R_{ij} which takes from one
    of the possible fitting basis (evolution, NN31IC, FLAVOUR) to the NN31IC basis,
    (sigma, g, v, v3, v8, t3, t8, cp), corresponding to the one used in NNPDF31.
    Denoting the rotation matrix as R_{ij} i is the flavour index and j is the evolution index.
    The evolution basis (NN31IC) is defined as
    cp = c + cbar = 2c
    and
    sigma = u + ubar + d + dbar + s + sbar + cp
    v = u - ubar + d - dbar + s - sbar + c - cbar
    v3 = u - ubar - d + dbar
    v8 = u - ubar + d - dbar - 2*s + 2*sbar
    t3 = u + ubar - d - dbar
    t8 = u + ubar + d + dbar - 2*s - 2*sbar

    If the input is already in the evolution basis it returns the identity.

    Parameters
    ----------
        flav_info: dict
            dictionary containing the information about each PDF (basis dictionary in the runcard)
        fitbasis: str
            name of the fitting basis

    Returns
    -------
        mat.transpose(): numpy matrix
            matrix performing the change of basis from fitbasis to NN31IC

    """
    if fitbasis == 'NN31IC':
        return np.identity(8)

    elif fitbasis == 'NN31PC':
        sng = {'sng': 1, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 'g': 0 }
        v =  {'sng': 0, 'v': 1, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 'g': 0 }
        v3 = {'sng': 0, 'v': 0, 'v3': 1, 'v8': 0, 't3': 0, 't8': 0, 'g': 0 }
        v8 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 1, 't3': 0, 't8': 0, 'g': 0 }
        t3 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 1, 't8': 0, 'g': 0 }
        t8 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 1, 'g': 0 }
        cp = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 'g': 0 }
        g = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 'g': 1 }

    elif fitbasis == 'FLAVOUR':
        sng = {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1, 'c': 2, 'g': 0 }
        v = {'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': 1, 'sbar': -1, 'c': 0, 'g': 0 }
        v3 = {'u': 1, 'ubar': -1, 'd': -1, 'dbar': 1, 's': 0, 'sbar': 0, 'c': 0, 'g': 0 }
        v8 = {'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1, 's': -2, 'sbar': 2, 'c': 0, 'g': 0 }
        t3 = {'u': 1, 'ubar': 1, 'd': -1, 'dbar': -1, 's': 0, 'sbar': 0, 'c': 0, 'g': 0 }
        t8 = {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': -2, 'sbar': -2, 'c': 0, 'g': 0 }
        cp = {'u': 0, 'ubar': 0, 'd': 0, 'dbar': 0, 's': 0, 'sbar': 0, 'c': 2, 'g': 0 }
        g = {'u': 0, 'ubar': 0, 'd': 0, 'dbar': 0, 's': 0, 'sbar': 0, 'c': 0, 'g': 1 }

    elif fitbasis == 'EVOL' or fitbasis == 'evolution':
        sng = {'sng': 1, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v = {'sng': 0, 'v': 1, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v3 = {'sng': 0, 'v': 0, 'v3': 1, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v8 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 1, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        t3 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 1, 't8': 0, 't15': 0, 'g': 0 }
        t8 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 1, 't15': 0, 'g': 0 }
        cp = {'sng': 0.25, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': -0.25, 'g': 0 }
        g = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 1 }

    elif fitbasis == 'PDF4LHC20':
        sng = {'sng': 1, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v = {'sng': 0, 'v': 1, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v3 = {'sng': 0, 'v': 0, 'v3': 1, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        v8 = {'sng': 0, 'v': 1, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        t3 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 1, 't8': 0, 't15': 0, 'g': 0 }
        t8 = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 1, 't15': 0, 'g': 0 }
        cp = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 0 }
        g = {'sng': 0, 'v': 0, 'v3': 0, 'v8': 0, 't3': 0, 't8': 0, 't15': 0, 'g': 1 }

    flist = [sng, g, v, v3, v8, t3, t8, cp]

    evol_basis = False
    mat = []
    for f in flist:
        for flav_dict in flav_info:
            flav_name = flav_dict["fl"]
            mat.append(f[flav_name])

    nflavs = len(flav_info)
    mat = np.asarray(mat).reshape(8, nflavs)
    # Return the transpose of the matrix, to have the first index referring to flavour
    return mat.transpose()
