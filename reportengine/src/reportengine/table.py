# -*- coding: utf-8 -*-
"""
Save generated tables.

Apply the decorators to functions returning a dataframe
and they will be saved in the correct path. If the providers return instead a
tuple (handle, df), the hanndle will be appended at the end of the file name
(but before the extension).

The tables are saved tab separated (i.e. the only reasonable way to write csv).

Otherwise, this works the same as figure.py.

Created on Mon May 16 11:36:13 2016

@author: Zahari Kassabov
"""
import logging
import re

#TODO: Get rid of explicit pandas dependency here?
import pandas as pd

from reportengine.formattingtools import spec_to_nice_name
from reportengine.utils import add_highlight
from reportengine.floatformatting import format_number

__all__ = ('table', 'tablegen')

log = logging.getLogger(__name__)

#http://pandas.pydata.org/pandas-docs/stable/internals.html#subclassing-pandas-data-structures
class Table(pd.DataFrame):

    _metadata = ['description', 'path']

    @classmethod
    def fromdf(cls, df, *, description=None, path=None):
        res = cls(df, copy=False)
        res.description = description
        res.path = path
        return res
    @property
    def _constructor(self):
        return Table

    @property
    def as_markdown(self):
        if self.size > 1e5:
            p = self.path.relative_to(self.path.parents[1])

            return f"> Very long [table]({p})"
        with pd.option_context('display.max_colwidth', -1):
            res = self.to_html(escape = False, border=0,
                               float_format=format_number)
        #See https://github.com/jgm/pandoc/issues/1841
        res = re.sub('\n\s+', '\n', res)
        return res



def prepare_path(*, spec, namespace,environment, **kwargs):
    name = spec_to_nice_name(namespace, spec)
    path = environment.table_folder / (name + '.csv')
    return {'path': path}

def savetable(df, path):
    """Final action to save figures, with a nice filename"""
    log.debug("Writing table %s" % path)
    df.to_csv(str(path), sep='\t', na_rep='nan')
    return Table.fromdf(df, path=path)

def savetablelist(dfs, path):
    """Final action to save lists of figures. It adds a numerical index as
    a suffix, for each figure in the generator."""

    res = []

    for i, df in enumerate(dfs):
        if isinstance(df, tuple):
            df, suffix = df
        else:
            suffix = str(i)
        tb_path = path.with_name(path.stem + suffix, path.suffix)
        tb = savetable(df, tb_path)
        res.append(tb)
    return res



@add_highlight
def table(f):
    """Save the resulting table as a tab separated csv file after
    it is generated"""
    f.prepare = prepare_path
    f.final_action = savetable
    return f
@add_highlight
def tablegen(f):
    """Save each table of the generator. See ``table``."""
    f.prepare = prepare_path
    f.final_action = savetablelist
    return f
