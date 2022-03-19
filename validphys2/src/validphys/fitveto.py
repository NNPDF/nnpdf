"""
fitveto.py

Module for the determination of passing fit replicas.

Current active vetoes:
   Positivity - Replicas with FitInfo.is_positive == False
   ChiSquared - Replicas with ChiSquared > nsigma_discard_chi2*StandardDev + Average
   ArclengthX - Replicas with ArcLengthX > nsigma_discard_arclength*StandardDev + Average
   Integrability - Replicas with IntegrabilityNumbers < integ_threshold
"""

import json
import logging
import numpy as np

log = logging.getLogger(__name__)

# Default thresholds for distribution vetos in units of standard deivations
NSIGMA_DISCARD_ARCLENGTH = 4.0
NSIGMA_DISCARD_CHI2 = 4.0
INTEG_THRESHOLD = 0.5


def distribution_veto(dist, prior_mask, nsigma_threshold):
    """ For a given distribution (a list of floats), returns a boolean mask
    specifying the passing elements. The result is a new mask of the elements that
    satisfy:

    value <=  mean + nsigma_threshold*standard_deviation

    Only points passing the prior_mask are
    considered in the average or standard deviation."""
    if sum(prior_mask) <= 1:
        return prior_mask
    dist = np.asarray(dist)
    passing = dist[prior_mask]
    average_pass = np.mean(passing)
    stderr_pass = np.std(passing)
    # NOTE that this has always not been abs
    # i.e replicas that are lower than the average by more than 4std pass
    return (dist - average_pass) <= nsigma_threshold * stderr_pass


def integrability_veto(dist, integ_threshold):
    """ For a given distribution (a list of floats), returns a boolean mask
    specifying the passing elements. The result is a new mask of the elements that
    satisfy:
    value <=  integ_threshold
    """
    dist = np.asarray(dist)
    return dist <= integ_threshold


def determine_vetoes(fitinfos: list, nsigma_discard_chi2: float, nsigma_discard_arclength: float, integ_threshold: float):
    """ Assesses whether replica fitinfo passes standard NNPDF vetoes
    Returns a dictionary of vetoes and their passing boolean masks.
    Included in the dictionary is a 'Total' veto.
    """

    # Setup distributions to veto upon: Make a dictionary {name: (values, threshold)}, where
    # values and threshold are to be filtered recusively as per ``distribution_veto``.
    # TODO ensure that all replicas have the same amount of arclengths
    distributions = {"ChiSquared": ([i.chi2 for i in fitinfos], nsigma_discard_chi2)}
    for i in range(0, len(fitinfos[0].arclengths)):
        distributions["ArcLength_" + str(i)] = (
            [j.arclengths[i] for j in fitinfos],
            nsigma_discard_arclength,
        )

    # Positivity veto
    posmask = np.array([replica.is_positive for replica in fitinfos], dtype=bool)
    # alphasmask = np.array([replica.alphas > 0.11666 for replica in fitinfos], dtype=bool)
    vetoes = {"Positivity": posmask}
    # vetoes["alphas"] = alphasmask
    # posmask = posmask & alphasmask
    total_mask = posmask.copy()

    # Integrability veto
    if len(fitinfos[0].integnumbers) == 0:
        log.warning(f"No integrability numbers in the fitinfo file")
    else:
        for i in range(0, len(fitinfos[0].integnumbers)):
            values = [j.integnumbers[i] for j in fitinfos] 
            key = "IntegNumber_" + str(i)
            vetoes[key] = integrability_veto(
                values, integ_threshold=integ_threshold)
    
    # Distribution vetoes
    while True:
        for key in distributions:
            values, threshold = distributions[key]
            vetoes[key] = distribution_veto(
                values, total_mask, nsigma_threshold=threshold
            )
        new_total_mask = np.all(list(vetoes.values()), axis=0)
        if sum(new_total_mask) == sum(total_mask):
            break
        total_mask = new_total_mask

    pass_chi2 = np.asarray(distributions["ChiSquared"][0])[total_mask]
    log.info(f"Passing average chi2: {np.mean(pass_chi2)}")

    vetoes["Total"] = total_mask
    return vetoes


def save_vetoes_info(veto_dict: dict, chi2_threshold, arclength_threshold, integ_threshold, filepath):
    """ Saves to file the chi2 and arclength thresholds used by postfit as well as veto
    dictionaries which contain information on which replicas pass each veto."""
    if filepath.exists():
        log.warning(f"Veto file {filepath} already exists. Overwriting file")
    with open(filepath, "w") as f:
        thresholds_dict = {
            "chi2_threshold": chi2_threshold,
            "arclength_threshold": arclength_threshold,
            "integrability_threshold": integ_threshold
        }
        veto_dict_tolist = {key: val.tolist() for key, val in veto_dict.items()}
        combined_dict = {**thresholds_dict, **veto_dict_tolist}
        json.dump(combined_dict, f)
