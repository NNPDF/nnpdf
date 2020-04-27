"""
Test to ensure the validphys.pseudodata.get_pseudodata action
correctly obtains the appropriate pseudodata for an n3fit fit.

To this end, a 10 replica fit has been uploaded named
`pseudodata_test_fit` obtained using 100 epochs, theoryID 162 and a
subset of DIS datasets. When this fit was performed, the `all_exp_infos`
was pickled and stored in `exp_infos.pickle` which is the benchmark
we use to ensure the action is working appropriately.
"""
import pickle
from importlib.resources import read_binary

import numpy as np
import pytest

from validphys.api import API
import validphys.tests.regressions

from reportengine.compat import yaml

EXAMPLE_RUNCARD = """fit: pseudodata_test_fit
pdf: pseudodata_test_fit

experiments:
  from_: fit

theory:
  from_: fit

t0pdfset:
  from_: datacuts

datacuts:
  from_: fit

theoryid:
  from_: theory

use_cuts: fromfit
"""


def test_pseudodata():
    exp_infos_bytes = read_binary(validphys.tests.regressions, "test_exp_infos.pickle")
    ns = yaml.safe_load(EXAMPLE_RUNCARD)
    # This is what all the fitted replicas saw
    exp_infos = pickle.loads(exp_infos_bytes)

    # We now need to convert these to postfit replicas
    fitted_indices = API.fitted_replica_indexes(**ns)
    fit_postfit_mapping = dict(enumerate(exp_infos, 1))
    exp_infos = [fit_postfit_mapping[i] for i in fitted_indices]

    pseudodata_info = API.get_pseudodata(**ns)

    # Loop over replicas
    for i, j in zip(exp_infos, pseudodata_info):
        # For each replica, loop over experiments
        for exp1, exp2 in zip(i, j):
            assert np.allclose(exp1["expdata"], exp2["expdata"])
            assert np.allclose(exp1["expdata_vl"], exp2["expdata_vl"])
