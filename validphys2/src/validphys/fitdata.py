# -*- coding: utf-8 -*-
"""
Utilities for loading data from fit folders
"""
import logging
from collections import namedtuple
from io import StringIO
import pathlib


import numpy as np
import yaml
from reportengine import collect
from reportengine.table import table
from reportengine.checks import make_argcheck, CheckError

from validphys.core import PDF
from validphys import checks
from validphys.plotoptions import get_info
from validphys import pdfgrids

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

ReplicaSpec = namedtuple('ReplicaSpec', ('index', 'path', 'info'))

FitInfo = namedtuple("FitInfo", ("nite", 'training', 'validation', 'chi2', 'pos_status', 'arclenghts'))
def load_fitinfo(replica_path, prefix):
    """Process the data in the ".fitinfo" file of a single replica."""
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

#TODO: Produce a more informative .sumrules file.
def load_sumrules(replica_path, prefix):
    """Load the values of the sum rules defined in
    ``validphys.pdfgrids.SUM_RULES`` from a given replica."""
    return np.loadtxt(replica_path/f'{prefix}.sumrules')[:len(pdfgrids.SUM_RULES)]

@checks.check_has_fitted_replicas
def replica_paths(fit):
    """Return the paths of all the replicas"""
    #Total number of members = number of replicas + 1
    l = len(PDF(fit.name))
    return [fit.path / 'nnfit' / f'replica_{index}' for index in range(1, l)]

def replica_data(fit, replica_paths):
    """Load the data from the fitinfo file of each of the replicas.
    The corresponding PDF set must be installed in the LHAPDF path.

    The included information is:

    ('nite', 'training', 'validation', 'chi2', 'pos_status', 'arclenghts')"""
    return [load_fitinfo(path, fit.name) for path in replica_paths]

def fit_sum_rules(fit, replica_paths):
    """Return a SumRulesGrid object with the sumrules for each replica as
    calculated by nnfit at the initial scale. This is the same object as
    the one produced by
    ``validphys.pdfgrids.sum_rules`` which is instead obtained from LHAPDF at
    a given energy"""
    res = np.zeros((len(pdfgrids.SUM_RULES),len(replica_paths)))
    for i, p in enumerate(replica_paths):
        res[:, i] = load_sumrules(p, fit.name)
    return pdfgrids.SumRulesGrid(*res)

@table
def fit_sum_rules_table(fit_sum_rules):
    return pdfgrids.sum_rules_table(fit_sum_rules)


fits_replica_data = collect('replica_data', ('fits',))

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
    """Given exactly two fits, print the datasets that are included in one "
    "but not in the other. If `print_common` is True, also print the datasets
    that are common."""
    m = match_datasets_by_name
    first,second = fits
    res = StringIO()
    if m.common and print_common:
        res.write("The following datasets are included in both `%s` and `%s`:\n\n" % (first, second))
        for k,v in m.common.items():
            info = get_info(v[0].commondata)
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    if m.first_only:
        res.write("The following datasets are included in `%s` but not in `%s`:\n\n"% (first,second))
        for k,v in m.first_only.items():
            info = get_info(v.commondata)
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    if m.second_only:
        res.write("The following datasets are included in `%s` but not in `%s`:\n\n"% (second,first))
        for k,v in m.second_only.items():
            info = get_info(v.commondata)
            res.write(' - %s\n' % info.dataset_label)
        res.write('\n')
    return res.getvalue()

print_dataset_differences.highlight = 'markdown'

@_assert_two_fits
def test_for_same_cuts(fits, match_datasets_by_name):
    """Given two fits, return a list of tuples `(first, second)`
    where `first` and `second` are
    DatasetSpecs that correspond to the same dataset but have different cuts,
    such that `first` is included in the first fit and `second` in the second.
    """
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

def print_different_cuts(fits, test_for_same_cuts):
    """Print a summary of the datasets that are included in both fits but have
    different cuts."""
    res = StringIO()
    first_fit, second_fit = fits
    if test_for_same_cuts:
        res.write("The following datasets are both included but have different kinematical cuts:\n\n")
        for (first, second) in test_for_same_cuts:
            info = get_info(first.commondata)
            total_points = len(first.commondata.load())
            res.write(" - %s:\n" % info.dataset_label)
            first_len = len(first.cuts.load()) if first.cuts else total_points
            second_len = len(second.cuts.load()) if second.cuts else total_points
            res.write("    * %s includes %d out of %d points.\n" % (first_fit, first_len, total_points))
            res.write("    * %s includes %d out of %d points.\n" % (second_fit, second_len, total_points))
        res.write('\n')


    return res.getvalue()

def _get_fitted_index(pdf, i):
    """Return the nnfit index for the replcia i"""
    p = pathlib.Path(pdf.infopath).with_name(f'{pdf.name}_{i:04d}.dat')
    with open(p) as f:
        it = yaml.load_all(f)
        metadata = next(it)
    return metadata['FromMCReplica']

@make_argcheck
def _check_has_replica_tags(pdf):
    """Check that the PDF has fitted index tags."""
    try:
        _get_fitted_index(pdf,1)
    except KeyError as e:
        raise CheckError("PDF replica file don't contain "
                         "the fitted replica tag.") from e

@_check_has_replica_tags
def fitted_replica_indexes(pdf):
    """Return nnfit index of replicas 1 to N."""
    return [_get_fitted_index(pdf,i) for i in range(1, len(pdf))]
