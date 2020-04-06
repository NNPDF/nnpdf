""" This module implements tools for computing convolutions. It reimplements
the relevant C++ functionality in pure Python (using numpy and pandas).

Note that currently no effort has been made to optimize these operations.
"""


import pandas as pd
import numpy as np

from validphys.pdfbases import evolution

__all__ = ('hadron_predictions', 'dis_predictions')

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


def hadron_predictions(pdf, loaded_fk):
    """Low level function to compute predictions from an hadronic
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
    """Low level function to compute predictions from a DIS
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
        >>> from validphys.convolution import dis_predictions
        >>> from validphys.fkparser import load_fktable
        >>> l = Loader()
        >>> pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
        >>> ds = l.check_dataset("HERAF2CHARM", theoryid=53)
        >>> table = load_fktable(ds.fkspecs[0]).with_cuts(ds.cuts)
        >>> from validphys.convolution import dis_predictions
        >>> dis_predictions(pdf, table).T.describe()
        data           15          16          17          18  ...          48          49          50          51
        count  100.000000  100.000000  100.000000  100.000000  ...  100.000000  100.000000  100.000000  100.000000
        mean     0.292386    0.269370    0.235357    0.198888  ...    0.099491    0.174955    0.086756    0.060422
        std      0.018171    0.015950    0.013248    0.011431  ...    0.006643    0.007640    0.005944    0.004877
        min      0.224384    0.214071    0.200633    0.167839  ...    0.075529    0.157506    0.066015    0.046424
        25%      0.284614    0.261566    0.228153    0.192455  ...    0.096242    0.170292    0.083439    0.057636
        50%      0.291761    0.267144    0.232952    0.196367  ...    0.100099    0.174754    0.087629    0.060917
        75%      0.301445    0.277181    0.242640    0.205627  ...    0.103017    0.179115    0.089806    0.062776
        max      0.356283    0.325746    0.281639    0.236154  ...    0.116030    0.196459    0.102557    0.074180
    """

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
