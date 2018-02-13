import numpy as np
import pytest
import itertools
from hypothesis import given
from hypothesis.strategies import (floats, integers, tuples, lists,
    booleans)
from hypothesis.extra.numpy import arrays,  array_shapes

from validphys.fitveto import (distribution_veto, determine_vetoes, NSIGMA_DISCARD)
from validphys.fitdata import FitInfo

shape1d = array_shapes(max_dims=1, min_side=1, max_side=1000)
nicefloats = floats(allow_nan=False, allow_infinity=False)
fitinfos = tuples(integers(min_value=1),
        nicefloats, nicefloats, nicefloats,
        booleans(),
        arrays(float, shape=7, elements=nicefloats)).map(FitInfo._make)


#Ignore over- and underflow warnings.
@pytest.mark.filterwarnings('ignore')
@given(arrays(float, shape=shape1d, elements=nicefloats))
def test_distribution_veto(arr):
    veto = distribution_veto(arr)
    masked = arr[veto]
    assert np.all(masked - np.mean(masked) <= NSIGMA_DISCARD*np.std(masked))
    assert sum(distribution_veto(masked)) == len(masked)


#The case where the list is empty is handled in postfit
@pytest.mark.filterwarnings('ignore')
@given(lists(fitinfos, min_size=1))
def test_determine_vetoes(fitinfos):
    vetoes = determine_vetoes(fitinfos)
    assert np.all(vetoes['Positivity'] == np.array([info.is_positive for info in fitinfos]))
    tot = vetoes['Total']
    assert all(np.all(tot & val == tot) for val in vetoes.values())
    # distribution_vetoes applied a second time should veto nothing
    if sum(tot) > 0:
        passing_fitinfos = list(itertools.compress(fitinfos, tot))
        second_vetoes = determine_vetoes(passing_fitinfos)
        assert sum(vetoes["Total"]) == sum(second_vetoes["Total"])
