# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
from reportengine.checks import make_argcheck, check

log = logging.getLogger(__name__)


@make_argcheck
def check_combocuts(combocuts: str):
    """Check combocuts content"""
    check(combocuts == 'NNPDF31',
          "Invalid combocut. Must be NNPDF31 (or implement it yourself).")


def make_dataset_dir(path):
    """check if results folder exists"""
    if path.exists():
        log.warning(f"Dataset output folder exists: {path} Overwritting contents")
    else:
        path.mkdir(exist_ok=True)


@check_combocuts
def filter(experiments, combocuts, theoryid, filter_path, t0set=None):
    """Apply filters to all datasets"""

    # Load experiments
    for experiment in experiments:
        for dataset in experiment.datasets:
            log.info('Filtering %s' % dataset)
            make_dataset_dir(filter_path / dataset.name)
            uncut_dataset = dataset.load()
            for idat in range(uncut_dataset.GetNData()):
                passKinCuts(uncut_dataset, idat)


def passKinCuts(dataset, idat):
    """Applies cuts as in C++ for NNPDF3.1 combo cuts"""
    pass
