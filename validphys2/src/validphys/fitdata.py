# -*- coding: utf-8 -*-
"""
Utilities for loading data from fit folders
"""
import pathlib
import logging
from collections import namedtuple

#import blessings
import numpy as np

from validphys.core import PDF
from validphys import checks
from validphys import lhaindex

#TODO: Add more stuff here as needed for postfit


#t = blessings.Terminal()
log = logging.getLogger(__name__)

#TODO: These should be the actual filenames, without the redundant prefix
REPLICA_FILES = (
'dat',
'fitinfo',
'params',
'preproc',
'sumrules',
)

LITERAL_FILES = (
'chi2exps.log',
'GAMin.log',
'nnfit.yml',
)




def check_results_path(path):
    path = pathlib.Path(path)
    assert path.is_dir(), 'Path is not a directory %s' % path
    assert (path / 'nnfit').is_dir(), 'Path "nnfit" is not a folder not in path'

ReplicaSpec = namedtuple('ReplicaSpec', ('index', 'path', 'info'))

FitInfo = namedtuple("FitInfo", ("nite", 'training', 'validation', 'chi2', 'pos_status', 'arclenghts'))
def load_fitinfo(replica_path, prefix):
    p = replica_path / (prefix + '.fitinfo')
    with p.open() as f:
        line = next(f)
        props = iter(line.split())
        nite = int(next(props))
        validation = float(next(props))
        training = float(next(props))
        chi2 = float(next(props))
        pos_status = next(props)
        line = next(f)
        arclenghts = np.fromstring(line, sep=' ')

    return FitInfo(nite, training, validation, chi2, pos_status, arclenghts)

@checks.check_has_fitted_replicas
def replica_data(fit):
    """Load the data from the fitinfo file of each of the replicas.
    The corresponding PDF set must be installed in the LHAPDF path.

    The included information is:

    ('nite', 'training', 'validation', 'chi2', 'pos_status', 'arclenghts')"""
    nreplicas = len(PDF(fit.name)) - 1
    infos = []
    for index in range(1, nreplicas+1):
        path = fit.path / 'nnfit' / ('replica_%d' % index)
        infos.append(load_fitinfo(path, fit.name))
    return infos