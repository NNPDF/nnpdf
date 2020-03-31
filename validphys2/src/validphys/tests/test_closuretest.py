"""
test_closuretest.py

contains some unit tests for closure test estimators

"""

import numpy as np

from validphys.closuretest import bias_dataset, variance_dataset

class TestResult:
    """class for testing base level estimators which expect a results object"""
    def __init__(self, central_value, rawdata=None):
        self.central_value = central_value
        self._rawdata = rawdata
        self.ndata = len(central_value)
    def __len__(self,):
        return self.ndata

N_DATA = 5
N_REPLICAS = 10

#TODO: make these fixtures?
ones_results = (None, TestResult(np.ones(N_DATA), np.ones((N_DATA, N_REPLICAS))))
twos_results = (None, TestResult(2*np.ones(N_DATA), 2*np.ones((N_DATA, N_REPLICAS))))

replicas = np.arange(N_REPLICAS)[np.newaxis, :]*np.ones((N_DATA, 1))
replicas_result = (None, TestResult(replicas.mean(axis=1), replicas))

def test_bias_function():
    bias_ones = bias_dataset(
        ones_results,
        [ones_results], # need list of length one to emulate collect
        None,
        None,
        np.identity(5)
    )
    assert np.allclose(0, bias_ones.bias)
    bias_one_two = bias_dataset(
        ones_results,
        [twos_results],
        None,
        None,
        np.identity(5)
    )
    assert np.allclose(N_DATA, bias_one_two.bias)

def test_variance_function():
    vardata = variance_dataset(
        ones_results,
        None,
        None,
        np.identity(5)
    )
    assert np.allclose(0, vardata.variance)
    var_reps = variance_dataset(
        replicas_result,
        None,
        None,
        np.identity(5)
    )
    # calc explicitly what variance should be
    expected = np.sum(((np.arange(N_REPLICAS) - 4.5)**2)*N_DATA/N_REPLICAS)
    assert np.allclose(expected, var_reps.variance)
