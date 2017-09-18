"""
#tableloader.py

Load from file some of the tables that validphys produces.
Contrary to `validphys.loader` this module consists of functions that take
absolute paths, and return mostly dataframes.
"""
import functools
import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

sane_load = functools.partial(pd.DataFrame.from_csv, sep='\t')

def fixup_header(df, head_index, dtype):
    """Set the type of the column index in place"""
    oldcols = df.columns
    good = oldcols.levels[head_index].map(dtype)
    newcols = oldcols.set_levels([*oldcols.levels[:head_index],
                                  good,
                                  *oldcols.levels[head_index+1:]])
    df.columns = newcols

def parse_exp_mat(filename):
    df = sane_load(filename, header=[0,1,2], index_col=[0,1,2])
    fixup_header(df, 2, int)
    return df

load_experiments_covmat = parse_exp_mat
load_experiments_invcovmat = parse_exp_mat

def load_perreplica_chi2_table(filename):
    """Load the output of ``perreplica_chi2_table``."""
    df = sane_load(filename, header=[0,1])
    fixup_header(df, 1, int)
    return df


def load_fits_computed_psedorreplicas_chi2(filename):
    """Load the output of ``fits_computed_psedorreplicas_chi2``"""
    return sane_load(filename, index_col=[0,1,2,3], header=[0,1,])


def load_fits_chi2_table(filename):
    return sane_load(filename, header=[0,1], index_col=[0,1])

def load_adapted_fits_chi2_table(filename):
    """Load the fits_chi2_table and adapt it in the way that suits the
    ``paramfits`` module."""
    df = load_fits_chi2_table(filename)
    f = lambda x: x[x.columns[0]]*x[x.columns[1]]
    df = df.groupby(axis=1, level=0).apply(f)
    df.columns = pd.MultiIndex.from_product([list(df.columns), ['chi2']])
    return df


def set_actual_column_level0(df, new_levels):
    """Set the fitst level of the index to new_levels.
    Pandas makes it criminally difficult to do this properly."""
    cols = df.columns
    levels = np.asarray(cols.levels[0])
    levels[cols.labels[0]] = new_levels
    cols.set_levels(levels, inplace=True, level=0)




#TODO: Find a better place for this function
def combine_pseudorreplica_tables(dfs, combined_names, blacklist_datasets=None):
    """Return a table in the same format as perreplica_chi2_table with th   e
    minimum value of the chi² for each batch of fits."""

    for df in dfs:
        set_actual_column_level0(df, combined_names)


    if blacklist_datasets:
        m = np.ones(df.shape[0], dtype=bool)
        for it in blacklist_datasets:
            dsmask = (df.index.get_level_values(1) != it)
            m &= dsmask
        if m.all():
            log.warning("Did not blacklist any dataset from the list {blacklisted_datasets}")
        else:
            df = df.loc[m]

    together = pd.concat(dfs, axis=1, keys=range(len(dfs)))

    total = together.loc[(slice(None), 'Total'), :]

    total_chis =  total.groupby(level=3).sum()

    #Note, asarray is needed because it ignores NANs otherwise.
    argmin = lambda x: pd.Series(np.argmin(np.asarray(x), axis=1), index=x.index)

    best_replicas = total_chis.groupby(axis=1, level=1).apply(argmin)
    gb = together.groupby(axis=1, level=1)

    def inner_select(df, indexes):
        return df.iloc[:,indexes[df.name]]


    def select_best_replicas(df):
        indexes = best_replicas[df.name]
        return df.groupby(level=3).apply(inner_select, indexes=indexes)

    res = gb.apply(select_best_replicas)
    res.index = res.index.droplevel(0)
    res.sort_index(inplace=True)

    #TODO: Why in earth did I decide to keep this?!
    res.columns = pd.MultiIndex.from_product((res.columns, ['chi2']))



    return res
