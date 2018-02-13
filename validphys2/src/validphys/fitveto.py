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


def distribution_veto(dist, prior_mask):
    """ For a given distribution (a list of floats), returns a boolean mask
    specifying the passing elements. Only points passing the prior_mask are
    considered in the average or standard deviation."""
    if sum(prior_mask) <= 1: return prior_mask
    dist = np.asarray(dist)
    passing = dist[prior_mask]
    average_pass = np.mean(passing)
    stderr_pass  = np.std(passing)
    # NOTE that this has always not been abs
    # i.e replicas that are lower than the average by more than 4std pass
    return (dist - average_pass) <= NSIGMA_DISCARD*stderr_pass


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
    total_mask = np.asarray(vetoes["Positivity"])

    # Distribution vetoes
    while True:
        for key in distributions:
            vetoes[key] = distribution_veto(distributions[key], total_mask)
        new_total_mask = np.all(list(vetoes.values()), axis=0)
        if sum(new_total_mask) == sum(total_mask): break
        total_mask = new_total_mask

    vetoes["Total"] = total_mask
    return vetoes
