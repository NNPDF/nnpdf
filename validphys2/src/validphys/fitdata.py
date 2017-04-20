# -*- coding: utf-8 -*-
"""
Utilities for loading data from fit folders
"""
import pathlib
import logging
from collections import namedtuple
from io import StringIO

from reportengine import collect
from reportengine.checks import make_argcheck, CheckError
import numpy as np

from validphys.core import PDF
from validphys import checks
from validphys import lhaindex
from validphys.plotoptions import get_infos

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

#Do collect in two parts so we get a list for each fit instead of a single list
all_datasets = collect('dataset', ('experiments', 'experiment'))
fits_datasets = collect('all_datasets', ('fits', 'fitinputcontext',))

@make_argcheck
def _assert_two_fits(fits):
    """Check that there are exatly two fits"""
    if len(fits) != 2:
        raise CheckError("Exactly two fits are required")

DatasetComp = namedtuple('DatasetComp', ('common', 'first_only', 'second_only'))

@_assert_two_fits
def match_datasets_by_name(fits, fits_datasets):
    """Return a tuple with common, first_only and second_only.
    The elements of the tuple are mappings where the keys are dataset names
    and the values are the two datasets contained in each fit for common, and
    the corresponfing dataset inclucded only in the first fit and only in the
    second fit."""

    firstds, secondds = [{ds.name: ds for ds in datasets} for datasets in fits_datasets]
    common_keys = firstds.keys() & secondds.keys()
    first_keys = firstds.keys() - secondds.keys()
    seccond_keys = secondds.keys() - firstds.keys()

    common = {k: (firstds[k], secondds[k]) for k in common_keys}
    first_only = {k: firstds[k] for k in first_keys}
    second_only = {k: secondds[k] for k in seccond_keys}
    return DatasetComp(common, first_only, second_only)


#TODO: Do we do md output here or that's for the templates?
def print_dataset_differences(fits, match_datasets_by_name,
                              print_common:bool=True):
    m = match_datasets_by_name
    first,second = fits
    res = StringIO()
    if m.common and print_common:
        res.write("The following datasets are included in both `%s` and `%s`:\n\n" % (first, second))
        for k,v in m.common.items():
            info = get_infos(v[0].commondata)[0]
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    if m.first_only:
        res.write("The following datasets were included in `%s` but not in `%s`:\n\n"% (first,second))
        for k,v in m.first_only.items():
            info = get_infos(v.commondata)[0]
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    if m.second_only:
        res.write("The following datasets were included in `%s` but not in `%s`:\n\n"% (second,first))
        for k,v in m.second_only.items():
            info = get_infos(v.commondata)[0]
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    return res.getvalue()


def test_for_same_cuts(fits, match_datasets_by_name):
    common = match_datasets_by_name.common
    first_fit, second_fit = fits
    res = []
    for ds, (first, second) in common.items():
        if first.cuts:
            c1 = first.cuts.load()
        else:
            c1 = ['Allpass']
        if second.cuts:
            c2 = second.cuts.load()
        else:
            c2 = ['Allpass']
        if not np.array_equal(c1, c2):
            msg = "Cuts for %s are not the same:\n%s:\n%s\n\n%s:\n%s" % (ds, first_fit, c1, second_fit, c2)
            log.info(msg)
            res.append((first,  second))
    return res

