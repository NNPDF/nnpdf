import numpy as np
import pytest
import itertools
from hypothesis import given
from hypothesis.strategies import (floats, integers, tuples, lists,
    booleans)
from hypothesis.extra.numpy import arrays,  array_shapes

from validphys.fitveto import distribution_veto, determine_vetoes
from validphys.fitveto import NSIGMA_DISCARD_ARCLENGTH, NSIGMA_DISCARD_CHI2
from validphys.fitdata import FitInfo

shape1d = array_shapes(max_dims=1, min_side=1, max_side=1000)
nicefloats = floats(allow_nan=False, allow_infinity=False)
fitinfos = tuples(integers(min_value=1),
        nicefloats, nicefloats, nicefloats,
        booleans(),
        arrays(float, shape=7, elements=nicefloats)).map(FitInfo._make)


thresholds = floats(min_value=1, max_value=10)
distributions = arrays(float, shape=shape1d, elements=nicefloats)
# Ignore over- and underflow warnings.
@pytest.mark.filterwarnings("ignore")
@given(distributions, thresholds)
def test_distribution_veto(arr, threshold):
    veto = distribution_veto(arr, np.ones_like(arr, dtype=bool), threshold)
    masked = arr[veto]
    assert np.all(masked - np.mean(arr) <= threshold * np.std(arr))


#The case where the list is empty is handled in postfit
@pytest.mark.filterwarnings('ignore')
@given(lists(fitinfos, min_size=1))
def test_determine_vetoes(fitinfos):
    vetoes = determine_vetoes(fitinfos, NSIGMA_DISCARD_CHI2, NSIGMA_DISCARD_ARCLENGTH)
    assert np.all(vetoes['Positivity'] == np.array([info.is_positive for info in fitinfos]))
    tot = vetoes['Total']
    assert all(np.all(tot & val == tot) for val in vetoes.values())
    single_replica_veto = determine_vetoes([fitinfos[0]], NSIGMA_DISCARD_CHI2, NSIGMA_DISCARD_ARCLENGTH)
    assert single_replica_veto['Total'][0] == single_replica_veto['Positivity'][0]
    # distribution_vetoes applied a second time should veto nothing
    if sum(tot) > 0:
        passing_fitinfos = list(itertools.compress(fitinfos, tot))
        second_vetoes = determine_vetoes(passing_fitinfos, NSIGMA_DISCARD_CHI2, NSIGMA_DISCARD_ARCLENGTH)
        assert sum(vetoes["Total"]) == sum(second_vetoes["Total"])
