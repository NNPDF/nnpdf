#!/usr/bin/env python

"""
This module contains a checks provider to be used by n3fit apps
"""

import hashlib
import logging
from pathlib import Path

import n3fit.checks

log = logging.getLogger(__name__)
MD5FK_FILENAME = "md5fk"


@n3fit.checks.check_consistent_basis
@n3fit.checks.wrapper_check_NN
@n3fit.checks.wrapper_hyperopt
@n3fit.checks.check_deprecated_options
@n3fit.checks.check_consistent_parallel
def n3fit_checks_action(
    *,
    genrep,
    data,
    theoryid,
    basis,
    fitbasis,
    sum_rules=True,
    parameters,
    save=None,
    load=None,
    load_weights_from_fit=None,
    load_weights_dict,
    hyperscan_config=None,
    hyperopt=None,
    kfold=None,
    tensorboard=None,
    parallel_models=True,
    double_precision=False,
    trial_specs=None,
    replicas=None,
):
    return


@n3fit.checks.check_eko_exists
def evolven3fit_checks_action(theoryid):
    return


def fktable_hasher(data, output_path):
    """Writes a hash of the FK-tables to a log file.
    This hash can be used to ensure whether two
    (supposedly identical) FK-tables of the same
    theory and dataset are numerically identical.
    """
    md5fk_path = output_path / MD5FK_FILENAME
    # Open a file to write in
    with open(md5fk_path, "w") as f:
        # Loop through the dataspecs object
        for dataset in data.datasets:
            fkspecs = dataset.fkspecs
            for fk in fkspecs:
                # Get the FK-table path
                fkpaths = fk.fkpath
                # For older theories, ensure we always pass a tuple
                if not isinstance(fkpaths, tuple):
                    fkpaths = (fk.fkpath,)
                # And then write the hash of the raw bytes
                for fkpath in fkpaths:
                    fkhash = hashlib.md5(fkpath.read_bytes()).hexdigest()
                    f.write(f"{Path(fkpath.stem).stem} {fkhash}\n")
        log.info(f"FK-table md5 hashes stored in {output_path}/md5fk")
