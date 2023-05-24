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

# NOTE:Considering the first columns as index by default (the index_col=0)
# is not particularly sane, but  turns out that it is advantageous for backward
# compatibility with the older DataFrame.from_csv method, that was employed
# previously.
sane_load = functools.partial(pd.read_csv, sep='\t', index_col=0)


class TableLoaderError(Exception):
    """Errors in the tableloader module."""

    pass


def fixup_header(df, head_index, dtype):
    """Set the type of the column index in place"""
    oldcols = df.columns
    good = oldcols.levels[head_index].map(dtype)
    newcols = oldcols.set_levels(
        [*oldcols.levels[:head_index], good, *oldcols.levels[head_index + 1 :]]
    )
    df.columns = newcols


def parse_data_cv(filename):
    """Useful for reading DataFrames with just one column."""
    df = sane_load(filename, index_col=[0, 1, 2])
    return df


def parse_exp_mat(filename):
    """Parse a dump of a matrix like experiments_covmat."""
    df = sane_load(filename, header=[0, 1, 2], index_col=[0, 1, 2])
    fixup_header(df, 2, int)
    return df


load_experiments_covmat = parse_exp_mat
load_experiments_invcovmat = parse_exp_mat


def load_perreplica_chi2_table(filename):
    """Load the output of ``perreplica_chi2_table``."""
    df = sane_load(filename, index_col=0, header=[0, 1])
    fixup_header(df, 1, int)
    return df


def load_fits_computed_pseudoreplicas_chi2(filename):
    """Load the output of ``fits_computed_psedorreplicas_chi2``"""
    return sane_load(filename, index_col=[0, 1, 2, 3], header=[0, 1])


def load_fits_chi2_table(filename):
    """Load the result of fits_chi2_tavle or similar."""
    return sane_load(filename, header=[0, 1], index_col=[0, 1])


def load_adapted_fits_chi2_table(filename):
    """Load the fits_chi2_table and adapt it in the way that suits the
    ``paramfits`` module. That is, return a table with the total chi² and
    another with the number of points."""
    df = load_fits_chi2_table(filename)
    ndatalabel = df.columns[0][1]
    dns = df.sort_index(axis=1).loc[:, pd.IndexSlice[:, ndatalabel]]
    if not (dns.apply(pd.Series.nunique, axis=1) == 1).all():
        raise TableLoaderError("Expecting all entries to have the same ndata")

    ndatas = dns.iloc[:, 0]

    f = lambda x: x[x.columns[0]] * x[x.columns[1]]
    df = df.groupby(axis=1, level=0).apply(f)
    df.columns = pd.MultiIndex.from_product([list(df.columns), ['chi2']])

    return ndatas, df


def set_actual_column_level0(df, new_levels):
    """Set the first level of the index to new_levels. Note:
    This is a separate function mostly because it breaks
    in every patch update of pandas."""
    cols = df.columns
    cols.set_levels(new_levels, inplace=True, level=0)


# TODO: Find a better place for this function
def combine_pseudoreplica_tables(
    dfs, combined_names, *, blacklist_datasets=None, min_points_required=2
):
    """Return a table in the same format as perreplica_chi2_table with th   e
    minimum value of the chi² for each batch of fits."""

    for df in dfs:
        set_actual_column_level0(df, combined_names)
        if blacklist_datasets:
            m = np.ones(df.shape[0], dtype=bool)
            for it in blacklist_datasets:
                dsmask = df.index.get_level_values(1) != it
                m &= dsmask
            if m.all():
                log.warning(f"Did not blacklist any dataset from the list {blacklist_datasets}")
            else:
                df = df.loc[m]

    together = pd.concat(dfs, axis=1, keys=range(len(dfs)))

    total = together.loc[(slice(None), 'Total'), :]

    total_chis = total.groupby(level=3).sum(min_count=1)

    def fixup_min_points(df):
        m = (~df.isnull()).sum(axis=1, min_count=1) >= min_points_required
        df[df[m].isnull()] = np.inf
        return df

    # The idea is: Set to inf the nans of the valid curves, so that we select
    # the minimum (which is not infinite).  Leave the bad nans as nans, so we
    # write nan always for those.
    total_chis = total_chis.groupby(axis=1, level=1).apply(fixup_min_points)

    # Note, asarray is needed because it ignores NANs otherwise.
    argmin = lambda x: pd.Series(np.argmin(np.asarray(x), axis=1), index=x.index)

    best_replicas = total_chis.groupby(axis=1, level=1).apply(argmin)
    gb = together.groupby(axis=1, level=1)

    def inner_select(df, indexes):
        return df.iloc[:, indexes[df.name]]

    def select_best_replicas(df):
        indexes = best_replicas[df.name]
        return df.groupby(level=3).apply(inner_select, indexes=indexes)

    res = gb.apply(select_best_replicas)
    res.index = res.index.droplevel(0)
    res.sort_index(inplace=True)

    # TODO: Why in earth did I decide to keep this?!
    res.columns = pd.MultiIndex.from_product((res.columns, ['chi2']))

    return res


def get_extrasum_slice(df, components):
    """Extract a slice of a table that has the components in the format that
    extra_sums expects."""
    df = pd.DataFrame(df)
    df.sort_index(inplace=True)
    total_token = ' Total'
    keys = [
        (c[: -len(total_token)], 'Total') if c.endswith(total_token) else (slice(None), c)
        for c in components
    ]
    locs = [flat for key in keys for flat in df.index.get_locs(key)]
    return df.iloc[locs, :]


# Define aliases for functions with spelling mistakes in their names which have now been corrected
# Do this so that old runcards still work
load_fits_computed_psedorreplicas_chi2 = load_fits_computed_pseudoreplicas_chi2
combine_pseudorreplica_tables = combine_pseudoreplica_tables
