""" This module implements tools for computing convolutions. It reimplements
the relevant C++ functionality in pure Python (using numpy and pandas).

The high level :py:func:`predictions` function can be used to extact theory
predictions for experimentally measured quantities, in a way that is directly
comparable with the C++ code::

    import numpy as np
    from validphys.api import API
    from validphys.convolution import predictions


    inp = {
        'fit': '181023-001-sc',
        'use_cuts': 'internal',
        'theoryid': 162,
        'pdf': 'NNPDF31_nlo_as_0118',
        'experiments': {'from_': 'fit'}
    }

    tb = API.experiment_result_table(**inp)

    all_datasets = [ds for e in API.experiments(**inp) for ds in e.datasets]

    pdf = API.pdf(**inp)


    all_preds = [predictions(ds, pdf) for ds in all_datasets]

    for ds, pred in zip(all_datasets, all_preds):
        cpp = tb.loc[(slice(None),ds.name), :]
        assert np.allclose(pred.values, cpp.values[:, 2:], atol=1e-7, rtol=1e-4)




Some variants such as :py:func:`central_predictions` and
:py:func:`linear_predictions` are useful for more specialized tasks.

These functions work with :py:class:`validphys.core.DatasetSpec` objects,
allowing to account for information on COMPOUND predictions and cuts. A lower
level interface which operates with :py:class:`validphys.coredata.FKTableData`
objects is also available.

Note that currently no effort has been made to optimize these operations.
"""
import operator
import functools

import pandas as pd
import numpy as np

from validphys.pdfbases import evolution
from validphys.fkparser import load_fktable


FK_FLAVOURS = evolution.to_known_elements(
    [
        "photon",
        "singlet",
        "g",
        "V",
        "V3",
        "V8",
        "V15",
        "V24",
        "V35",
        "T3",
        "T8",
        "T15",
        "T24",
        "T35",
    ]
)

NFK = len(FK_FLAVOURS)


def _asy(a, b):
    return (a - b) / (a + b)


def _smn(a, b, c, d):
    return (a + b) / (c + d)

def _com(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t):
    return (a + b + c + d + e + f + g + h + i + j) / ( k + l + m + n + o + p + q + r + s + t)

def _smt(a, b, c, d, e, f, g, h, i, j):
    return (a + b + c + d + e + f + g + h + i + j)

def _id(a):
    return a


OP = {
    "RATIO": operator.truediv,
    "ASY": _asy,
    "ADD": operator.add,
    "SMN": _smn,
    "COM": _com,
    "SMT": _smt,
    "NULL": _id,
}

class PredictionsRequireCutsError(Exception): pass

def _predictions(dataset, pdf, fkfunc):
    """Combine data on all the FKTables in the database according to the
    reduction operation defined therein. Dispatch the kind of predictions (for
    all replicas, central, etc) according to the provided ``fkfunc``, which
    should have the same interface as e.g. ``fk_predictions``.
    """
    opfunc = OP[dataset.op]
    if dataset.cuts is None:
        raise PredictionsRequireCutsError(
            "FKTables do not always generate predictions for some datapoints "
            "which are usually cut. Loading predictions without cuts can "
            "therefore produce predictions whose shape doesn't match the uncut "
            "commondata and is not supported."
        )
    cuts = dataset.cuts.load()
    all_predictions = [
        fkfunc(load_fktable(fk).with_cuts(cuts), pdf) for fk in dataset.fkspecs
    ]
    return opfunc(*all_predictions)


def predictions(dataset, pdf):
    """"Compute theory predictions for a given PDF and dataset. Information
    regading the dataset, on cuts, CFactors and combinations of FKTables is
    taken into account to construct the predictions.

    The result should be comparable to experimental predictions implemented in
    CommonData.

    Parameters
    ----------
    dataset : validphys.core.DatasetSpec
        The dataset containing information on the partonic cross section.
    pdf : validphys.core.PDF
        The PDF set to use for the convolutions.

    Returns
    -------
    df : pandas.DataFrame
        A dataframe corresponding to the hadronic prediction for each data
        point for the PDF members. The index of the dataframe corresponds to
        the selected data points, based on the dataset :ref:`cuts <filters>`. The
        columns correspond to the selected PDF members in the LHAPDF set, which
        depend on the PDF error type (see
        :py:meth:`validphys.core.PDF.grid_values_index`)

    Examples
    --------
    Obtain descriptive statistics over PDF replicas for each of the three
    points in the ATLAS ttbar dataset:


    >>> from validphys.loader import Loader
    >>> l = Loader()
    >>> ds = l.check_dataset('ATLASTTBARTOT', theoryid=53)
    >>> from validphys.convolution import predictions
    >>> pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
    >>> preds = predictions(ds, pdf)
    >>> preds.T.describe()
    data            0           1           2
    count  100.000000  100.000000  100.000000
    mean   161.271292  231.500367  767.816844
    std      2.227304    2.883497    7.327617
    min    156.638526  225.283254  750.850250
    25%    159.652216  229.486793  762.773527
    50%    161.066965  231.281248  767.619249
    75%    162.620554  233.306836  772.390286
    max    168.390840  240.287549  786.549380


    """
    return _predictions(dataset, pdf, fk_predictions)


def central_predictions(dataset, pdf):
    """Same as :py:func:`predictions` but computing the predictions for the
    central member of the PDF set only. For Monte Carlo PDFs, this is a faster
    alternative to computing the central predictions as the average of the
    replica predictions (although a small approximation is involved in the case
    of hadronic predictions).
    """
    return _predictions(dataset, pdf, central_fk_predictions)


def linear_predictions(dataset, pdf):
    """Same as :py:func:`predictions` but computing *linearized* predictions.
    These are the same as ``predictions`` for DIS, but truncates to the terms
    that are linear in the difference between each member and the central
    value for hadronic predictions.

    This approximation is generally a very good approximation in that yields
    differences that are much smaller that the PDF uncertainty.
    """
    return _predictions(dataset, pdf, linear_fk_predictions)


def fk_predictions(loaded_fk, pdf):
    """Low level function to compute predictions from a
    FKTable.

    Parameters
    ----------
    loaded_fk : validphys.coredata.FKTableData
        The FKTable corresponding to the partonic cross section.
    pdf :  validphys.core.PDF
        The PDF set to use for the convolutions.

    Returns
    -------
    df : pandas.DataFrame
        A dataframe corresponding to the hadronic prediction for each data
        point for the PDF members. The index of the dataframe corresponds to
        the selected data points (use
        :py:meth:`validphys.coredata.FKTableData.with_cuts` to filter out
        points). The columns correspond to the selected PDF members in the
        LHAPDF set, which depend on the PDF error type (see
        :py:meth:`validphys.core.PDF.grid_values_index`)

    Notes
    -----
    This function operates on a single FKTable, while the prediction for an
    experimental quantity generally involves several. Use
    :py:func:`predictions` to compute those.

    Examples
    --------

        >>> from validphys.loader import Loader
        >>> from validphys.convolution import hadron_predictions
        >>> from validphys.fkparser import load_fktable
        >>> l = Loader()
        >>> pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
        >>> ds = l.check_dataset('ATLASTTBARTOT', theoryid=53, cfac=('QCD',))
        >>> table = load_fktable(ds.fkspecs[0])
        >>> hadron_predictions(table, pdf)
                     1           2           3           4    ...         97          98          99          100
        data                                                  ...
        0     176.688118  170.172930  172.460771  173.792321  ...  179.504636  172.343792  168.372508  169.927820
        1     252.682923  244.507916  247.840249  249.541798  ...  256.410844  247.805180  242.246438  244.415529
        2     828.076008  813.452551  824.581569  828.213508  ...  838.707211  826.056388  810.310109  816.824167

    """
    if loaded_fk.hadronic:
        return hadron_predictions(loaded_fk, pdf)
    else:
        return dis_predictions(loaded_fk, pdf)


def central_fk_predictions(loaded_fk, pdf):
    """Same as :py:func:`fk_predictions`, but computing predictions for the
    central PDF member only."""
    if loaded_fk.hadronic:
        return central_hadron_predictions(loaded_fk, pdf)
    else:
        return central_dis_predictions(loaded_fk, pdf)


def linear_fk_predictions(loaded_fk, pdf):
    """Same as :py:func:`predictions` for DIS, but compute linearized
    predictions for hadronic data, using :py:func:`linear_hadron_predictions`.
    """
    if loaded_fk.hadronic:
        return linear_hadron_predictions(loaded_fk, pdf)
    else:
        return dis_predictions(loaded_fk, pdf)


def _gv_hadron_predictions(loaded_fk, gv1func, gv2func=None):
    """Compute hadronic convolutions between the loaded FKTable
    and the PDF evaluation functions `gv1func` and `gv2func`.
    These must have the same interface as
    :py:meth:`validphys.pdfbases.evolution.grid_values`, but without the PDF
    argument.

    If gv2func is not given, then gv1func will be used for the second PDF,
    with the grid being evaluated only once.
    """
    xgrid = loaded_fk.xgrid
    Q = loaded_fk.Q0
    sigma = loaded_fk.sigma

    # Generate gid values for all flavours in the evolution basis, in the
    # expected order.
    #
    # Squeeze to remove the dimension over Q.
    gv1 = gv1func(qmat=[Q], vmat=FK_FLAVOURS, xmat=xgrid).squeeze(-1)
    if gv2func is not None:
        gv2 = gv2func(qmat=[Q], vmat=FK_FLAVOURS, xmat=xgrid).squeeze(-1)
    else:
        gv2 = gv1

    # The hadronic FK table columns are indexes into the NFK*NFK table of
    # possible flavour combinations of the two PDFs, with the convention of
    # looping first of the first index and the over the second: If the flavour
    # index of the first PDF is ``i`` and the second is ``j``, then the column
    # value in the FKTable is ``i*NFK + j``. This can easily be inverted using
    # the ``np.indices``, which is used here to map the column index to i and
    # j.
    fm = sigma.columns
    all_fl_indices_1, all_fl_indices_2 = np.indices((NFK, NFK))
    # The flavour indices of the first and second PDF for each combination
    # (column) are the columns indexing into the flattened indices.
    fl1 = all_fl_indices_1.ravel()[fm]
    fl2 = all_fl_indices_2.ravel()[fm]
    # Once we have the flavours, shape the PDF grids as appropriate for the
    # convolution below: We are left with two tensor of shape
    # ``nmembers * len(sigma.columns) * nx`` such that the pairs of flavours of the two
    # combinations correspond to the combination encoded in the FKTable.
    expanded_gv1 = gv1[:, fl1, :]
    expanded_gv2 = gv2[:, fl2, :]
    # Create a luminosity tensor holding the value f1(x1)*f2(x2) for all
    # possible x1-x2 combinations (f1, f2, x1, x2)
    luminosity = np.einsum("ijk, ijl->ijkl", expanded_gv1, expanded_gv2)

    def appl(df):
        # x1 and x2 are encoded as the first and second index levels.
        xx1 = df.index.get_level_values(1)
        xx2 = df.index.get_level_values(2)
        # take the active combinations from the luminosity tensor
        partial_lumi = luminosity[..., xx1, xx2]
        return pd.Series(np.einsum("ijk,kj->i",partial_lumi, df.values))

    return sigma.groupby(level=0).apply(appl)


def _gv_dis_predictions(loaded_fk, gvfunc):
    xgrid = loaded_fk.xgrid
    Q = loaded_fk.Q0
    sigma = loaded_fk.sigma
    # The column indexes are indices into the FK_FLAVOURS list above.
    fm = sigma.columns
    # Squeeze to remove the dimension over Q.
    gv = gvfunc(qmat=[Q], vmat=FK_FLAVOURS[fm], xmat=xgrid).squeeze(-1)

    def appl(df):
        # x is encoded as the first index level.
        xind = df.index.get_level_values(1)
        return pd.Series(np.einsum("ijk,kj->i", gv[:, :, xind], df.values))

    return sigma.groupby(level=0).apply(appl)


def hadron_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`fk_predictions` for hadronic observables."""
    gv = functools.partial(evolution.grid_values, pdf=pdf)
    res = _gv_hadron_predictions(loaded_fk, gv)
    res.columns = pdf.grid_values_index
    return res


def central_hadron_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`central_fk_predictions` for hadronic
    observables."""
    gv = functools.partial(evolution.central_grid_values, pdf=pdf)
    return _gv_hadron_predictions(loaded_fk, gv)


def linear_hadron_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`linear_fk_predictions` for hadronic
    observables. Specifically this computes:

        central_value ⊗ FK ⊗ (2 * replica_values - central_value)

    which is the linear expansion of the hadronic observable in the difference
    between each replica and the central value, ``replica_values -
    central_value``
    """
    gv1 = functools.partial(evolution.central_grid_values, pdf=pdf)

    def gv2(*args, **kwargs):
        replica_values = evolution.grid_values(pdf, *args, **kwargs)
        central_value = evolution.central_grid_values(pdf, *args, **kwargs)
        return 2 * replica_values - central_value

    res = _gv_hadron_predictions(loaded_fk, gv1, gv2)
    res.columns = pdf.grid_values_index
    return res


def dis_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`fk_predictions` for DIS observables."""
    gv = functools.partial(evolution.grid_values, pdf=pdf)
    res = _gv_dis_predictions(loaded_fk, gv)
    res.columns = pdf.grid_values_index
    return res


def central_dis_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`central_fk_predictions` for DIS
    observables."""
    gv = functools.partial(evolution.central_grid_values, pdf=pdf)
    return _gv_dis_predictions(loaded_fk, gv)


def _positivity_predictions(posdataset, pdf, fkfunc):
    """Implentation of :py:func:`_predictions` but for positivity
    datasets."""
    return fkfunc(load_fktable(posdataset.fkspec), pdf)


def positivity_predictions(posdataset, pdf):
    """Implementation of :py:func:`predictions` but for positivity
    datasets."""
    return _positivity_predictions(posdataset, pdf, fk_predictions)


def linear_positivity_predictions(posdataset, pdf):
    """Implmentation of :py:func:`linear_predictions` but for positivity
    datasets."""
    return _positivity_predictions(posdataset, pdf, linear_fk_predictions)


def central_positivity_predictions(posdataset, pdf):
    """Implementation as :py:func:`central_predictions`, but for postivity datasets
    """
    return _positivity_predictions(posdataset, pdf, central_fk_predictions)
