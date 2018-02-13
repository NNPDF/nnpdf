"""
fitveto.py

Module for the determination of passing fit replicas.

Current active vetoes:
   Positivity - Replicas with FitInfo.is_positive == False
   ChiSquared - Replicas with ChiSquared > NSIGMA_DISCARD*StandardDev + Average
   ArclengthX - Replicas with ArcLengthX > NSIGMA_DISCARD*StandardDev + Average
"""

import numpy as np

# Threshold for distribution vetos
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


def determine_vetoes(fitinfos: list):
    """ Assesses whether replica fitinfo passes standard NNPDF vetoes
    Returns a dictionary of vetoes and their passing boolean masks.
    Included in the dictionary is a 'Total' veto.
    """
    # Setup distributions to veto upon
    # TODO ensure that all replicas have the same amount of arclengths
    distributions = {"ChiSquared": [i.chi2 for i in fitinfos]}
    for i in range(0, len(fitinfos[0].arclengths)):
        distributions["ArcLength_"+str(i)] = [j.arclengths[i] for j in fitinfos]

    # Positivity veto
    vetoes = {"Positivity": [replica.is_positive for replica in fitinfos]}

    # Distribution vetoes
    for key in distributions:
        vetoes[key] = distribution_veto(distributions[key])

    vetoes["Total"] = np.ones(len(fitinfos), dtype=bool)
    for key in vetoes:
        vetoes["Total"] = np.asarray([x & y for (x, y) in
            zip(vetoes["Total"], vetoes[key])])
    return vetoes
