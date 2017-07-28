"""
#tableloader.py

Load from file some of the tables that validphys produces.
Contrary to `validphys.loader` this module consists of functions that take
absolute paths, and return mostly dataframes.
"""
import functools

import pandas as pd

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
    df = sane_load(filename, header=[0,1])
    fixup_header(df, 1, int)
    return df


def load_fits_computed_psedorreplicas_chi2(filename):
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

#TODO: Find a better place for this function
def combine_pseudorreplica_tables(dfs, combined_names):
    """Return a table in the same format as perreplica_chi2_table with the
    minimum value of the chi² for each batch of fits."""
    for df in dfs:
        df.columns = df.columns.set_levels(combined_names, level=0)



    together = pd.concat(dfs, axis=1, keys=range(len(dfs)))
    res = together.min(axis=1, level=1, skipna=False)
    res.columns = dfs[0].columns


    #Kind of horrible. What we do is to select the total, collapse the replica
    #index (min could really be anything) and sum across the index.
    total = together.loc[(slice(None), 'Total'), :]
    #Note this is NOT the same as the sum of the minima.
    total_chis =  total.groupby(level=3).sum().min(axis=1, level=1, skipna=False)
    ntotal = total.min(level=(0,1,2)).index.get_level_values(2).get_values().sum()
    total_index = pd.MultiIndex.from_product(
            [['Total'], ['Total'], [ntotal], total_chis.index],
            names = res.index.names)
    total_chis.index = total_index
    total_chis.columns = res.columns
    res = pd.concat([res, total_chis])



    return res
