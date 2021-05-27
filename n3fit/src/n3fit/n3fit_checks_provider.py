#!/usr/bin/env python

"""
This module contains a checks provider to be used by n3fit apps
"""

import n3fit.checks

@n3fit.checks.can_run_parallel_replicas
@n3fit.checks.check_consistent_basis
@n3fit.checks.wrapper_check_NN
@n3fit.checks.wrapper_hyperopt
@n3fit.checks.check_deprecated_options
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
    parallel_models=1,
):
    # can_run_parallel_replicas(genrep, parameters, hyperopt, parallel_models)
    # check_consistent_basis(sum_rules, fitbasis, basis, theoryid)
    # wrapper_check_NN(basis, tensorboard, save, load, parameters)
    # wrapper_hyperopt(hyperopt, hyperscan, genrep, data)
    # check_deprecated_options(fitting)
    pass
