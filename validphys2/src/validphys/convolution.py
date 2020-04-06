""" This module implements tools for computing convolutions. It reimplements
the relevant C++ functionality in pure Python (using numpy and pandas).

Note that currently no effort has been made to optimize these operations.
"""


import pandas as pd
import numpy as np

from validphys.pdfbases import evolution

__all__ = ('predictions',)

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


def predictions(pdf, loaded_fk):
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

    gv = evolution.grid_values(pdf=pdf, qmat=[Q], vmat=FK_FLAVOURS, xmat=xgrid).squeeze(
        -1
    )
    fm = sigma.columns
    fl1, fl2 = np.indices((NFK, NFK))
    fl1 = fl1.ravel()[fm]
    fl2 = fl2.ravel()[fm]
    gv1fl = gv[:, fl1, :]
    gv2fl = gv[:, fl2, :]

    def appl(df):
        xx1 = df.index.get_level_values(1)
        xx2 = df.index.get_level_values(2)
        gv1 = gv1fl[:, :, xx1]
        gv2 = gv2fl[:, :, xx2]
        return pd.Series(np.einsum("ijk,ijk,kj->i", gv1, gv2, df.values), index=index)

    return sigma.groupby(level=0).apply(appl)


def dis_predictions(pdf, loaded_fk):
    xgrid = loaded_fk.xgrid
    Q = loaded_fk.Q0
    sigma = loaded_fk.sigma
    index = pdf.grid_values_index
    fm = sigma.columns
    gv = evolution.grid_values(
        pdf=pdf, qmat=[Q], vmat=FK_FLAVOURS[fm], xmat=xgrid
    ).squeeze(-1)

    def appl(df):
        xind = df.index.get_level_values(1)
        return pd.Series(np.einsum("ijk,kj->i", gv[:, :, xind], df.values), index=index)

    return sigma.groupby(level=0).apply(appl)
