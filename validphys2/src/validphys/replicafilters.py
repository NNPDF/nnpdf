# -*- coding: utf-8 -*-
"""
Filters for postfit.
"""
import logging

import numpy as np

from validphys.utils import split_by
from validphys.fitdata import LITERAL_FILES, REPLICA_FILES

log = logging.getLogger(__name__)

NSIGMA_DISCARD = 4



def valid_replica(path, prefix):
    if not path.is_dir():
        return False
    existing_files = set(path.iterdir())
    valid = (all(path/f in existing_files for f in LITERAL_FILES) and
            all(path/(prefix+'.'+f) for f in REPLICA_FILES))

    if not valid:
        log.warn("Found invalid replica %s" % path)
    return valid

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


