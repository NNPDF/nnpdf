#!/usr/bin/env python

"""
This module contains a checks provider to be used by n3fit apps
"""

import n3fit.checks


@n3fit.checks.check_consistent_basis
@n3fit.checks.fktable_hasher
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
    output_path,
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
@n3fit.checks.fktable_hasher
def fktable_hasher(data, output_path):
    return
