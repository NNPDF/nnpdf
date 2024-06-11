"""
test_closuretest.py

contains some unit tests for closure test module

"""

import numpy as np
from unittest.mock import Mock

from validphys.api import API
from validphys.tests.conftest import SINGLE_DATASET, THEORYID, PDF
from validphys.closuretest import bias_dataset, variance_dataset
from validphys.closuretest.multiclosure import (
    internal_multiclosure_dataset_loader,
    internal_multiclosure_data_loader,
)
from validphys.results import ThPredictionsResult


class TestResult:
    """class for testing base level estimators which expect a results object"""

    def __init__(self, central_value, rawdata=None):
        self.central_value = central_value
        self.rawdata = rawdata
        self.error_members = rawdata
        self.ndata = len(central_value)
        self.sqrtcovmat = np.identity(self.ndata)

    def __len__(self):
        return self.ndata


N_DATA = 5
N_REPLICAS = 10

# TODO: make these fixtures?
# these are proxies for results tuples of data and theory
ones_results = 2 * [TestResult(np.ones(N_DATA), np.ones((N_DATA, N_REPLICAS)))]
twos_results = 2 * [TestResult(2 * np.ones(N_DATA), 2 * np.ones((N_DATA, N_REPLICAS)))]

replicas = np.arange(N_REPLICAS)[np.newaxis, :] * np.ones((N_DATA, 1))
replicas_result = 2 * [TestResult(replicas.mean(axis=1), replicas)]

DATASET_INPUT = {"dataset_input": SINGLE_DATASET, "theoryid": THEORYID, "use_cuts": "internal"}

DATASET_INPUTS = {"dataset_inputs": [SINGLE_DATASET], "theoryid": THEORYID, "use_cuts": "internal"}

MOCK_FIT = Mock()
MOCK_FIT.name = "Mock fit spec"
MOCK_FIT.path = "mock/fit/path"
MOCK_FIT.label = MOCK_FIT.name

MOCK_FITS = [MOCK_FIT]


def test_internal_multiclosure_dataset_loader():
    """
    Test the internal_multiclosure_dataset_loader
    """
    dataset = API.dataset(**DATASET_INPUT)

    t0_covmat_from_systematics = API.t0_covmat_from_systematics(
        **{**DATASET_INPUT, "use_t0": True, "t0pdfset": PDF}
    )

    fits_pdf = [API.pdf(pdf=PDF)]
    multiclosure_underlyinglaw = API.pdf(pdf=PDF)
    fits = MOCK_FITS

    loader = internal_multiclosure_dataset_loader(
        dataset, fits_pdf, multiclosure_underlyinglaw, fits, t0_covmat_from_systematics
    )

    assert loader.closures_th is not None
    assert type(loader.closures_th) == list
    assert loader.law_th is not None
    assert type(loader.law_th) == ThPredictionsResult
    assert loader.covmat is not None
    assert type(loader.covmat) == np.ndarray
    assert loader.sqrt_covmat is not None
    assert type(loader.sqrt_covmat) == np.ndarray


def test_internal_multiclosure_data_loader():
    """
    Test the internal_multiclosure_data_loader
    """
    data = API.data(**DATASET_INPUTS)

    t0_covmat_from_systematics = API.t0_covmat_from_systematics(
        **{**DATASET_INPUT, "use_t0": True, "t0pdfset": PDF}
    )

    fits_pdf = [API.pdf(pdf=PDF)]
    multiclosure_underlyinglaw = API.pdf(pdf=PDF)
    fits = MOCK_FITS

    loader = internal_multiclosure_data_loader(
        data, fits_pdf, multiclosure_underlyinglaw, fits, t0_covmat_from_systematics
    )

    assert loader.closures_th is not None
    assert type(loader.closures_th) == list
    assert loader.law_th is not None
    assert type(loader.law_th) == ThPredictionsResult
    assert loader.covmat is not None
    assert type(loader.covmat) == np.ndarray
    assert loader.sqrt_covmat is not None
    assert type(loader.sqrt_covmat) == np.ndarray


def test_bias_function():
    bias_ones = bias_dataset(
        ones_results, [ones_results], None, None  # need list of length one to emulate collect
    )
    assert np.allclose(0, bias_ones.bias)
    bias_one_two = bias_dataset(ones_results, [twos_results], None, None)
    assert np.allclose(N_DATA, bias_one_two.bias)


def test_variance_function():
    vardata = variance_dataset(ones_results, None, None)
    assert np.allclose(0, vardata.variance)
    var_reps = variance_dataset(replicas_result, None, None)
    # calc explicitly what variance should be
    expected = np.sum(((np.arange(N_REPLICAS) - 4.5) ** 2) * N_DATA / N_REPLICAS)
    assert np.allclose(expected, var_reps.variance)
