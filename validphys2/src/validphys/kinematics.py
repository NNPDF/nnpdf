# -*- coding: utf-8 -*-
"""
Provides information on the kinematics involved in the data.

Uses the PLOTTING file specification.
"""
import logging

import pandas as pd
import numpy as np

from reportengine import collect
from reportengine.table import table

from validphys import plotoptions

log = logging.getLogger(__name__)

#TODO: Take commondata and cuts instead of dataset
def kinlimits(dataset):
    """Return a mapping with 'dataset', 'var', 'min', 'max' where 'var'
    corresponds to the x variable as defined in the PLOTTING specification and
    'min' and 'max' are the minimum and maximum of this variable for all the
    datapoints.
    """
    infos = plotoptions.get_infos(dataset)
    if len(infos)>1:
        log.info("Reading the first info for dataset %s "
        "and ignoring the others.", dataset)
    info = infos[0]
    kintable = plotoptions.kitable(dataset, info)
    varlabel = info.xlabel
    if info.x == 'idat':
        lmin, lmax = np.nan, np.nan
    else:
        lmin, lmax = kintable[info.x].min(), kintable[info.x].max()
    return {'dataset': dataset, 'var': varlabel, 'min':lmin, 'max':lmax}

all_kinlimits = collect(kinlimits, ('experiments', 'experiment'))

@table
def all_kinlimits_table(all_kinlimits):
    """Return a table with the kinematic limits for the datasets in all
    the experiments."""
    return pd.DataFrame(all_kinlimits, columns=['dataset', 'var', 'min', 'max'])
