import numpy as np
import pandas as pd
from scipy.special import xlogy
from typing import Tuple
from tqdm import tqdm
import argparse

class Unweight:
    def __init__(self, weights: np.ndarray, *chis: Tuple[int, np.ndarray]) -> None:
        """
        Initialize the Unweight class.

        Args:
            weights (np.ndarray): Array of weights.
            *chis (Tuple[int, np.ndarray]): Variable number of tuples containing power `n` and `chi` array.

        Raises:
            AssertionError: If lengths of `chi` arrays are not consistent.
        """
        self.chis = chis
        
        length = len(chis[0][1])
        for chi in chis:
            assert len(chi[1]) == length, "Not all chis have the same length!"

        if weights is None:
            self.weights = np.ones(length)
        else:
            self.weights = weights

    def entropy(self, p1, p2):
        """
        Calculate the entropy between two probability distributions.

        Args:
            p1 (np.ndarray): Probability distribution 1.
            p2 (np.ndarray): Probability distribution 2.

        Returns:
            float: Entropy value.
        """
        log = xlogy(p1, p1/p2)
        log[np.isnan(log)] = 0
        entropy = np.sum(log)
        return entropy
    
    def reweight(self, i: int = 0, thresh: float = 1e-12) -> None:
        """
        Perform reweighting.

        Args:
            i (int): Index of the chi array to reweight.
            thresh (float, optional): Threshold value for setting small weights to zero. Defaults to 1e-12.
        """

        n, chi = self.chis[i]
        exp = (n-1)*np.log(chi) - 1/2*np.power(chi,2.0)
        self.reweights = np.exp(exp - np.mean(exp))
        self.reweights = len(self.reweights)*self.reweights/np.sum(self.reweights)
        self.reweights[self.reweights <= thresh] = 0
        self.reprobs = self.reweights/len(self.reweights)

    def unweight(self, Np: int) -> None:
        """
        Perform unweighting.

        Args:
            Np (int): Number of points.
        """
        pcum = np.zeros(len(self.reweights) + 1)
        pcum[1:] = np.cumsum(self.reprobs)
        unweights = np.zeros(len(self.reweights), dtype="int")
        for k in range(len(self.reweights)):
            for j in range(Np):
                condition_one = j/Np - pcum[k] >= 0
                condition_two = pcum[k+1] - j/Np >= 0
                if condition_one and condition_two:
                    unweights[k] += 1
        
        self.unweights = unweights
        self.unprobs = unweights/np.sum(unweights)

    def effective_replicas(self, weights: np.ndarray) -> int:
        """
        Calculate the effective number of replicas.

        Args:
            weights (np.ndarray): Array of weights.

        Returns:
            int: Effective number of replicas.
        """
        N = len(weights)
        Neff = int(np.exp(-1/N*np.sum(xlogy(weights,weights/N))))
        return Neff

    def optimize(self, thresh: float, earlystopping: bool = True):
        """
        Optimize the unweighting process based on entropy threshold.

        Args:
            thresh (float): Entropy threshold value.
            earlystopping (bool, optional): Whether to stop optimization early if threshold is reached. Defaults to True.

        Returns:
            Tuple[np.ndarray, np.ndarray, int]: Tuple containing arrays of Nps, entropies, and optimal Np value.
        """
        Nps = np.logspace(1, np.log10(len(self.weights))+1, 50, dtype=np.int64)
        entropies = np.zeros(len(Nps))
        for i in tqdm(range(len(Nps))):
            self.unweight(Nps[i])
            entropies[i] = self.entropy(self.unprobs, self.reprobs)
            if entropies[i] <= thresh and earlystopping:
                loc = i
                break

        if i == len(Nps)-1:
            try:
                loc = np.where(entropies <= thresh)[0][0]
            except:
                print("Failed minimisation procedure! Defaulting to lowest entropy.")
                loc = -1
                
        Nopt = Nps[loc]

        return Nps, entropies, Nopt
    
def main(chi2, N, store = True):
    u = Unweight(None, (N, np.sqrt(chi2)))
    u.reweight()
    Neff = u.effective_replicas(u.reweights)
    u.unweight(Neff)

    weights = pd.DataFrame()
    weights["unweight"] = u.unweights
    weights["reweight"] = u.reweights
    weights["nrep"] = np.arange(1, len(weights)+1)
    weights = weights.set_index("nrep", drop = True)

    if store:
        weights.to_csv("weights.csv")

    return weights

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unweighting using chi squared values.")
    parser.add_argument("chi2_name", help = "Add the filename of the chi2 dataset (.csv)")
    parser.add_argument("N", help = "Add the amount of experimental datapoints that the chi2 is based on")
    args = parser.parse_args()
    chi2 = pd.read_csv(args.chi2_name).values

    main(chi2, args.N)