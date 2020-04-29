""" This module implements tools for computing convolutions. It reimplements
the relevant C++ functionality in pure Python (using numpy and pandas).

The high level :py:func:`predictions` function can be used to extact theory
predictions for experimentally measured quantities, in a way that is directly
comparable with the C++ code::

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



Note that currently no effort has been made to optimize these operations.
"""
import operator

import pandas as pd
import numpy as np

from validphys.pdfbases import evolution
from validphys.fkparser import load_fktable

__all__ = ('fk_predictions', 'predictions', 'dis_predictions', 'hadron_predictions')


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


def _id(a):
    return a


OP = {
    "RATIO": operator.truediv,
    "ASY": _asy,
    "ADD": operator.add,
    "SMN": _smn,
    "NULL": _id,
}


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
    opfunc = OP[dataset.op]
    cuts = dataset.cuts.load() if dataset.cuts is not None else None
    all_predictions = [
        fk_predictions(load_fktable(fk).with_cuts(cuts), pdf) for fk in dataset.fkspecs
    ]
    return opfunc(*all_predictions)


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


def hadron_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`fk_predictions` for hadronic observables."""
    xgrid = loaded_fk.xgrid
    Q = loaded_fk.Q0
    sigma = loaded_fk.sigma
    index = pdf.grid_values_index

    # Generate gid values for all flavours in the evolution basis, in the
    # expected order.
    #
    # Squeeze to remove the dimension over Q.
    gv = evolution.grid_values(pdf=pdf, qmat=[Q], vmat=FK_FLAVOURS, xmat=xgrid).squeeze(
        -1
    )
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
    expanded_gv1 = gv[:, fl1, :]
    expanded_gv2 = gv[:, fl2, :]

    def appl(df):
        # x1 and x2 are encoded as the first and second index levels.
        xx1 = df.index.get_level_values(1)
        xx2 = df.index.get_level_values(2)
        gv1 = expanded_gv1[:, :, xx1]
        gv2 = expanded_gv2[:, :, xx2]
        return pd.Series(np.einsum("ijk,ijk,kj->i", gv1, gv2, df.values), index=index)

    return sigma.groupby(level=0).apply(appl)


def dis_predictions(loaded_fk, pdf):
    """Implementation of :py:func:`fk_predictions` for DIS observables."""
    xgrid = loaded_fk.xgrid
    Q = loaded_fk.Q0
    sigma = loaded_fk.sigma
    index = pdf.grid_values_index
    # The column indexes are indices into the FK_FLAVOURS list above.
    fm = sigma.columns
    # Squeeze to remove the dimension over Q.
    gv = evolution.grid_values(
        pdf=pdf, qmat=[Q], vmat=FK_FLAVOURS[fm], xmat=xgrid
    ).squeeze(-1)

    def appl(df):
        # x is encoded as the first index level.
        xind = df.index.get_level_values(1)
        return pd.Series(np.einsum("ijk,kj->i", gv[:, :, xind], df.values), index=index)

    return sigma.groupby(level=0).apply(appl)
