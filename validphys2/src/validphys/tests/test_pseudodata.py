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
from validphys.pseudodata import training_validation_pseudodata
import validphys.tests.regressions

from reportengine.checks import CheckError
from reportengine.compat import yaml
from reportengine.resourcebuilder import ResourceError

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


@pytest.fixture(
    scope="session",
    params=[1, pytest.param(None, marks=pytest.mark.linux)],
)
def setup_dicts(request):
    n_process_config = dict(NPROC=request.param)
    exp_infos_bytes = read_binary(validphys.tests.regressions, "test_exp_infos.pickle")
    ns = yaml.safe_load(EXAMPLE_RUNCARD)
    # This is what all the fitted replicas saw
    exp_infos = pickle.loads(exp_infos_bytes)

    # We now need to convert these to postfit replicas
    fitted_indices = API.fitted_replica_indexes(**ns)
    fit_postfit_mapping = dict(enumerate(exp_infos, 1))
    exp_infos = [fit_postfit_mapping[i] for i in fitted_indices]

    pseudodata_info = API.get_pseudodata(**ns, **n_process_config)

    return exp_infos, pseudodata_info


def test_read_fit_pseudodata():
    data_generator = API.read_fit_pseudodata(
      fit="NNPDF31_nnlo_as_0118_DISonly_pseudodata",
      use_cuts="fromfit"
    )

    # Only bother checking the first ten replicas
    for _ in range(10):
      data, tr_idx, val_idx = next(data_generator)
      # Check the training and validation index are disjoint
      assert set(tr_idx).isdisjoint(set(val_idx))
      # Check the union is equal to the full dataset
      assert all(tr_idx.union(val_idx) == data.index)

    with pytest.raises(ResourceError) as e_info:
        API.read_fit_pseudodata(
          fit="NNPDF31_nnlo_as_0118_DISonly",
          use_cuts="nocuts"
        )
        assert isinstance(e_info.__cause__, CheckError)


def test_pseudodata(setup_dicts):
    exp_infos, pseudodata_info = setup_dicts
    # Loop over replicas
    for i, j in zip(exp_infos, pseudodata_info):
        # For each replica, loop over experiments
        for exp1, exp2 in zip(i, j):
            assert np.allclose(exp1["expdata"], exp2["expdata"])
            assert np.allclose(exp1["expdata_vl"], exp2["expdata_vl"])


def test_pseudodata_generator(setup_dicts):
    exp_infos, pseudodata_info = setup_dicts
    gen = training_validation_pseudodata(pseudodata_info)
    for i, j in enumerate(gen):
        continue
    # There is only one postfit replica in this fit
    assert i == 0
    # The training and validation split should be disjoint
    assert set(j["trdata"].index).isdisjoint(j["vldata"].index)
