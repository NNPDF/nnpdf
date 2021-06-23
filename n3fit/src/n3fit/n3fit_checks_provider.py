#!/usr/bin/env python

"""
This module contains a checks provider to be used by n3fit apps
"""

import n3fit.checks

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
    hyperscan=None,
    hyperopt=None,
    tensorboard=None,
    parallel_models=False,
    same_trvl_per_replica=False
):
    return
