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

parse_experiments_covmat = parse_exp_mat
parse_experiments_invcovmat = parse_exp_mat

def parse_perreplica_chi2_table(filename):
    df = sane_load(filename, header=[0,1])
    fixup_header(df, 1, int)
    return df