# -*- coding: utf-8 -*-
"""
Provides information on the kinematics involved in the data.

Uses the PLOTTING file specification.
"""
import logging

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table

from validphys import plotoptions

log = logging.getLogger(__name__)

def kinlimits(commondata, cuts, use_cuts, use_kinoverride:bool=True):
    """Return a mapping conaining the number of fitted and used datapoints,
    as well as the label, minimum and maximum value for each of the three
    kinematics. If ``use_kinoverride`` is set to False, the PLOTTING files
    will be ignored and the kinematics will be interpred based on the process
    type only. If use_cuts is False, the information on the total number of
    points will be displayed, instead of the fitted ones."""
    infos = plotoptions.get_infos(commondata, cuts=None, use_plotfiles=use_kinoverride)
    if len(infos)>1:
        log.info("Reading the first info for dataset %s "
        "and ignoring the others.", commondata)
    info = infos[0]
    kintable = plotoptions.kitable(commondata, info)
    ndata = len(kintable)
    if cuts:
        kintable = kintable.ix[cuts.load()]
        nfitted = len(kintable)
    elif use_cuts:
        nfitted = len(kintable)
    else:
        nfitted = '-'

    d = {'dataset': commondata, '$N_{data}$':ndata, '$N_{fitted}$':nfitted}
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
