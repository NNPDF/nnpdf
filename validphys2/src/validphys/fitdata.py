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
from validphys.utils import split_by
from validphys import checks
from validphys import lhaindex



#t = blessings.Terminal()
log = logging.getLogger(__name__)

NSIGMA_DISCARD = 4

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

def valid_replica(path, prefix):
    if not path.is_dir():
        return False
    existing_files = set(path.iterdir())
    valid = (all(path/f in existing_files for f in LITERAL_FILES) and
            all(path/(prefix+'.'+f) for f in REPLICA_FILES))

    if not valid:
        log.warn("Found invalid replica %s" % path)
    return valid


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
        training = float(next(props))
        validation = float(next(props))
        chi2 = float(next(props))
        pos_status = next(props)
        line = next(f)
        arclenghts = np.fromstring(line, sep=' ')

    return FitInfo(nite, training, validation, chi2, pos_status, arclenghts)


#TODO: Make postfit know about which replicas it has selected
@checks.make_check
def check_has_fitted_replicas(ns, **kwargs):
    name, path = ns['fit']
    postfit_path = path/'nnfit'/'postfit.log'
    if not postfit_path.exists():
        raise checks.CheckError("Fit {name} does not appear to be completed. "
        "Expected to find file {postfit_path}".format(**locals()))
    pdf = PDF(name)
    if not pdf.isinstalled:
        raise checks.CheckError("The PDF corresponding to the fit, '%s', "
        "needs to be "
        "installed in LHAPDF (i.e. copied to %s)."%
        (name, lhaindex.get_lha_datapath()))

@check_has_fitted_replicas
def replica_data(fit):
    nreplicas = len(PDF(fit.name)) - 1
    infos = []
    for index in range(1, nreplicas+1):
        path = fit.path / 'nnfit' / ('replica_%d' % index)
        infos.append(load_fitinfo(path, fit.name))
    return infos



def filter_positivity(repspec):
    if repspec.info.pos_status == 'POS_PASS':
        return True
    else:
        log.debug("Replica %s does not pass the positivity veto" % repspec.path.name)
        return False


def filter_chi2(repspecs):
    chi2 = np.array([rep.info.chi2 for rep in repspecs])
    m, s = np.mean(chi2), np.std(chi2)
    mask = np.abs(chi2 - m) < NSIGMA_DISCARD*s

    log.info("Mean chi² is: %.2f" % m)



    good, bad = split_by(repspecs, mask)

    if bad:
        log.info("Discarding %d replicas because of bad chi²." % len(bad))
    else:
        log.info("All replicas pass the chi² veto")
    if log.isEnabledFor(logging.DEBUG):
        for spec in bad:
            log.debug("Removed replica %s due to bad chi2." % spec.path.name)

    return good, bad

def filter_arclength(repspecs):
    alens = np.array([spec.info.arclenghts for spec in repspecs])
    m, s = np.mean(alens, axis=0), np.std(alens, axis=0)
    #Keep replica (axis 0) if all flavours (axis 1)
    #are within the acceptable range
    mask = (np.abs(alens - m) < NSIGMA_DISCARD*s).all(axis=1)

    good, bad = split_by(repspecs, mask)

    if bad:
        log.info("Discarding %d replicas because of bad arclength." % len(bad))
    else:
        log.info("All replicas pass the chi² veto")
    if log.isEnabledFor(logging.DEBUG):
        for spec in bad:
            log.debug("Removed replica %s due to bad chi2." % spec.path.name)

    return good, bad

