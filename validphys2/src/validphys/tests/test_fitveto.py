import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import (floats, integers, tuples, lists,
    booleans)
from hypothesis.extra.numpy import arrays,  array_shapes

from validphys.fitdata import FitInfo


#TODO: Import these from some sensible location

NSIGMA_DISCARD = 4

def distribution_veto(dist):
    """ For a given distribution (a list of floats), returns a boolean mask
    specifying the passing elements """
    dist = np.asarray(dist)
    replica_mask = np.ones_like(dist, dtype=bool)
    while True:
        passing = dist[replica_mask]
        average_pass = np.mean(passing)
        stderr_pass  = np.std(passing)
        # NOTE that this has always not been abs
        # i.e replicas that are lower than the average by more than 4std pass
        new_mask = (dist - average_pass) < NSIGMA_DISCARD*stderr_pass
        if sum(new_mask) == sum(replica_mask): break
        replica_mask = new_mask
    return replica_mask

def determine_vetoes(fitinfo: list):
    """ Assesses whether replica fitinfo passes standard NNPDF vetoes
    Returns a dictionary of vetoes and their passing boolean masks.
    Included in the dictionary is a 'Total' veto.
    """
    # Setup distributions to veto upon
    # TODO ensure that all replicas have the same amount of arclengths
    distributions = {"ChiSquared": [i.chi2 for i in fitinfo]}
    for i in range(0, len(fitinfo[0].arclengths)):
        distributions["ArcLength_"+str(i)] = [j.arclengths[i] for j in fitinfo]

    # Positivity veto
    vetoes = {"Positivity": [replica.is_positive for replica in fitinfo]}

    # Distribution vetoes
    for key in distributions:
        vetoes[key] = distribution_veto(distributions[key])

    # Determine total veto
    vetoes["Total"] = np.ones(len(fitinfo), dtype=bool)
    for key in vetoes:
        vetoes["Total"] = np.asarray([x & y for (x, y) in
            zip(vetoes["Total"], vetoes[key])])
    return vetoes


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
    assert np.all(veto - np.mean(veto) <= NSIGMA_DISCARD*np.std(veto))


#TODO: Fix the case where the list is empty somewhere
@pytest.mark.filterwarnings('ignore')
@given(lists(fitinfos, min_size=1))
def test_determine_vetoes(fitinfos):
    vetoes = determine_vetoes(fitinfos)
    assert np.all(vetoes['Positivity'] == np.array([info.is_positive for info in
            fitinfos]))
    tot = vetoes['Total']
    assert all(np.all(tot & val == tot) for val in vetoes.values())

