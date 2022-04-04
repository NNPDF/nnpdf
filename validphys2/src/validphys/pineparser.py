"""
    Loader for the pineappl-based FKTables

    The FKTables for pineappl have ``pineappl.lz4`` and can be utilized
    directly with the ``pineappl`` cli as well as read with ``pineappl.fk_table``
"""
import numpy as np
import pandas as pd
from validphys.coredata import FKTableData

# TODO: move back pineko things here


def get_yaml_information(yaml_file, theorypath):
    from pineko.parser import get_yaml_information as pineko_yaml, pineappl_to_fktable
    grids_folder = theorypath / "pineappls"
    return pineko_yaml(yaml_file, grids_folder)


def pineappl_reader(fkspec):
    """
    Receives a fkspec, which contains the appropiate references to the pineapplgrid so that
    pineko can parse it.

    The output of pineko contains all necessary information to create a dataframe
    with a multiindex of (data, x1, x2) -or only (data, x) for DIS-
    and a column per luminosity channel
    """
    from pineko.parser import get_yaml_information as pineko_yaml, pineappl_to_fktable
    fkdata = pineappl_to_fktable(fkspec.metadata, fkspec.fkpath)
    # The result is hadronic if the dataframe contains (data, x1 and x2)
    hadronic = fkdata.hadronic
    xgrid = fkdata.xgrid

    # Now loop over the fktables inside this observable and make them into dataframes
    partial_fktables = []
    for raw_fktable, lumi_columns, data_idx in fkdata:
        lf = len(lumi_columns)

        # Create the multi-index for the dataframe
        xi = np.arange(len(xgrid))

        if hadronic:
            idx = pd.MultiIndex.from_product([data_idx, xi, xi], names=["data", "x1", "x2"])
        else:
            idx = pd.MultiIndex.from_product([data_idx, xi], names=["data", "x"])

        # Now concatenate (data, x1, x2) and move the flavours to the columns
        df_fktable = raw_fktable.swapaxes(0,1).reshape(lf, -1).T
        partial_fktables.append(pd.DataFrame(df_fktable, columns=lumi_columns, index=idx))

    # Finallly concatenate all fktables, sort by flavours and fill any holes
    sigma = pd.concat(partial_fktables, sort=True, copy=False).fillna(0.0)

    # the number of data points is whatever the index of the last row is
    # since the indices can be shifted and we could have holes in the middle
    ndata = sigma.index.get_level_values(0).unique()[-1] + 1

    return FKTableData(
        sigma=sigma,
        ndata=ndata,
        Q0=fkdata.q0,
        metadata=fkspec.metadata,
        hadronic=hadronic,
        xgrid=xgrid,
        protected=fkdata.protected,
    )
