# -*- coding: utf-8 -*-
"""
Provides information on the kinematics involved in the data.

Uses the PLOTTING file specification.
"""
from collections import namedtuple
import logging

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table
from reportengine.checks import check_positive

from validphys import plotoptions

log = logging.getLogger(__name__)



@check_positive('titlelevel')
def describe_kinematics(commondata, titlelevel:int=1):
    """Output a markdown text describing the stored metadata for a given
    commondata.

    titlelevel can be used to control the header level of the title.
    """
    import inspect
    cd = commondata
    info = plotoptions.get_info(cd)
    proc = cd.load().GetProc(0)
    src = inspect.getsource(info.kinematics_override.xq2map)
    titlespec = '#'*titlelevel
    return (f"""
{titlespec} {cd}

{info.dataset_label}

Stored data:

 - Process type: **{proc}** ({info.process_description})

 - variables:
     * k1: {info.kinlabels[0]}
     * k2: {info.kinlabels[1]}
     * k3: {info.kinlabels[2]}



Map:

```python
{src}
```

""")

describe_kinematics.highlight = 'markdown'


nfittedlabel = '$N_{fitted}$'
ndatalabel = '$N_{data}$'

def kinlimits(commondata, cuts, use_cuts, use_kinoverride:bool=True):
    """Return a mapping conaining the number of fitted and used datapoints,
    as well as the label, minimum and maximum value for each of the three
    kinematics. If ``use_kinoverride`` is set to False, the PLOTTING files
    will be ignored and the kinematics will be interpred based on the process
    type only. If use_cuts is False, the information on the total number of
    points will be displayed, instead of the fitted ones."""
    info = plotoptions.get_info(commondata, cuts=None, use_plotfiles=use_kinoverride)

    kintable = plotoptions.kitable(commondata, info)
    ndata = len(kintable)
    if cuts:
        kintable = kintable.ix[cuts.load()]
        nfitted = len(kintable)
    elif use_cuts:
        nfitted = len(kintable)
    else:
        nfitted = '-'

    d = {'dataset': commondata, ndatalabel:ndata, nfittedlabel:nfitted}
    for i, key in enumerate(['k1', 'k2', 'k3']):
        kmin = kintable[key].min()
        kmax = kintable[key].max()
        label = info.kinlabels[i]
        d[key] = label
        d[key + ' min'] = kmin
        d[key + ' max'] = kmax
    return d

all_kinlimits = collect(kinlimits, ('experiments', 'experiment'))

@table
def all_kinlimits_table(all_kinlimits, use_kinoverride:bool=True):
    """Return a table with the kinematic limits for the datasets in all
    the experiments. If the PLOTTING overrides are not used, the information on
    sqrt(k2) will be displayed."""

    table = pd.DataFrame(all_kinlimits,
        columns=['dataset', '$N_{data}$', '$N_{fitted}$',
        'k1', 'k1 min', 'k1 max', 'k2', 'k2 min', 'k2 max', 'k3', 'k3 min', 'k3 max'
    ])

    #We really want to see the square root of the scale
    if not use_kinoverride:
        table['k2'] = 'sqrt(' + table['k2'] + ')'
        table['k2 min'] = np.sqrt(table['k2 min'])
        table['k2 max'] = np.sqrt(table['k2 max'])
        #renaming the columns is overly complicated
        cols = list(table.columns)
        cols[6:9]  = ['sqrt(k2)', 'sqrt(k2) min', 'sqrt(k2) max']
        table.columns = cols


    return table

def total_fitted_points(all_kinlimits_table)->int:
    """Print the total number of fitted points in a given set of experiments"""
    tb = all_kinlimits_table
    return int(tb[nfittedlabel].sum())


XQ2Map = namedtuple('XQ2Map', ('experiment', 'commondata', 'fitted', 'masked'))

def xq2map_with_cuts(experiment, commondata, cuts):
    """Return two (x,Q²) tuples: one for the fitted data and one for the
    cut data. If `display_cuts` is false or all data passes the cuts, the second
    tuple will be empty."""
    info = plotoptions.get_info(commondata)
    kintable = plotoptions.kitable(commondata, info)
    if cuts:
        mask = cuts.load()
        boolmask = np.zeros(len(kintable), dtype=bool)
        boolmask[mask] = True
        fitted_kintable = kintable.ix[boolmask]
        masked_kitable = kintable.ix[~boolmask]
        xq2fitted =  plotoptions.get_xq2map(fitted_kintable, info)
        xq2masked = plotoptions.get_xq2map(masked_kitable, info)
        return XQ2Map(experiment, commondata, xq2fitted, xq2masked)
    fitted_kintable = plotoptions.get_xq2map(kintable, info)
    empty = (np.array([]), np.array([]))
    return XQ2Map(experiment, commondata, fitted_kintable, empty)

experiments_xq2map = collect(xq2map_with_cuts, ('experiments', 'experiment'))
