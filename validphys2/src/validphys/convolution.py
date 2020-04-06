""" This module implements tools for computing convolutions. It reimplements
the relevant C++ functionality in pure Python (using numpy and pandas).

Note that currently no effort has been made to optimize these operations.
"""
import operator

import pandas as pd
import numpy as np

from validphys.pdfbases import evolution
from validphys.fkparser import load_fktable

__all__ = ('fk_predictions', 'predictions')


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


def predictions(pdf, dataset):
    opfunc = OP[dataset.op]
    cuts = dataset.cuts.load() if dataset.cuts is not None else None
    all_predictions = [
        fk_predictions(pdf, load_fktable(fk).with_cuts(cuts)) for fk in dataset.fkspecs
    ]
    return opfunc(*all_predictions)


def fk_predictions(pdf, loaded_fk):
    """Low level function to compute predictions from an
    observable.

    Parameters
    ----------
    pdf :  validphys.core.PDF
        The PDF set to use for the convolutions.
    loaded_fk : validphys.coredata.FKTableData
        The FKTable corresponding to the partonic cross section.

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

    Examples
    --------

        >>> from validphys.loader import Loader
        >>> from validphys.convolution import hadron_predictions
        >>> from validphys.fkparser import load_fktable
        >>> l = Loader()
        >>> pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
        >>> ds = l.check_dataset('ATLASTTBARTOT', theoryid=53, cfac=('QCD',))
        >>> table = load_fktable(ds.fkspecs[0])
        >>> hadron_predictions(pdf ,table)
                     1           2           3           4    ...         97          98          99          100
        data                                                  ...
        0     176.688118  170.172930  172.460771  173.792321  ...  179.504636  172.343792  168.372508  169.927820
        1     252.682923  244.507916  247.840249  249.541798  ...  256.410844  247.805180  242.246438  244.415529
        2     828.076008  813.452551  824.581569  828.213508  ...  838.707211  826.056388  810.310109  816.824167

    """
    if loaded_fk.hadronic:
        return hadron_predictions(pdf, loaded_fk)
    else:
        return dis_predictions(pdf, loaded_fk)


def hadron_predictions(pdf, loaded_fk):
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


def dis_predictions(pdf, loaded_fk):
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
