import numpy as np
import pandas as pd
from scipy.special import xlogy
from typing import Tuple
import argparse
import matplotlib.pyplot as plt


class Unweight:
    def __init__(self, weights: np.ndarray, *chis: Tuple[int, np.ndarray]) -> None:
        """
        Initialize the Unweight class.

        Parameters
        ----------
        weights : np.ndarray
            Array of weights.
        *chis : Tuple[int, np.ndarray]
            Variable number of tuples containing power `n` and `chi` array.

        Raises
        ------
        AssertionError
            If lengths of `chi` arrays are not consistent.
        """
        self.chis = chis

        length = len(chis[0][1])
        for chi in chis:
            assert len(chi[1]) == length, "Not all chis have the same length!"

        if weights is None:
            self.weights = np.ones(length)
        else:
            self.weights = weights

    def entropy(self, p1: np.ndarray, p2: np.ndarray) -> float:
        """
        Calculate the entropy between two probability distributions.

        Parameters
        ----------
        p1 : np.ndarray
            Probability distribution 1.
        p2 : np.ndarray
            Probability distribution 2.

        Returns
        -------
        float
            Entropy value.
        """
        log = xlogy(p1, p1 / p2)
        log[np.isnan(log)] = 0
        entropy = np.sum(log)
        return entropy

    def reweight(self, i: int = 0, thresh: float = 1e-12) -> None:
        """
        Perform reweighting.

        Parameters
        ----------
        i : int, optional
            Index of the chi array to reweight.
        thresh : float, optional
            Threshold value for setting small weights to zero. Defaults to 1e-12.
        """

        n, chi = self.chis[i]
        exp = (n - 1) * np.log(chi) - 1 / 2 * np.power(chi, 2.0)
        self.reweights = np.exp(exp - np.mean(exp))
        self.reweights = len(self.reweights) * self.reweights / np.sum(self.reweights)
        self.reweights[self.reweights <= thresh] = 0
        self.reprobs = self.reweights / len(self.reweights)

    def unweight(self, Np: int) -> None:
        """
        Perform unweighting.

        Parameters
        ----------
        Np : int
            Number of points.
        """
        pcum = np.zeros(len(self.reweights) + 1)
        pcum[1:] = np.cumsum(self.reprobs)
        unweights = np.zeros(len(self.reweights), dtype="int")
        for k in range(len(self.reweights)):
            for j in range(Np):
                condition_one = j / Np - pcum[k] >= 0
                condition_two = pcum[k + 1] - j / Np >= 0
                if condition_one and condition_two:
                    unweights[k] += 1

        self.unweights = unweights
        self.unprobs = unweights / np.sum(unweights)

    def effective_replicas(self, weights: np.ndarray, thresh: float = 1e-12) -> int:
        """
        Calculate the effective number of replicas.

        Parameters
        ----------
        weights : np.ndarray
            Array of weights.
        thresh : float, optional
            Threshold value neglecting small weights. Defaults to 1e-12.

        Returns
        -------
        int
            Effective number of replicas.
        """
        N = len(weights)
        weights = weights[weights > thresh]
        Neff = int(np.exp(-1 / N * np.sum(xlogy(weights, weights / N))))
        return Neff

    def iterate(
        self, thresh: float = 0, earlystopping: bool = True
    ) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Itarate the unweighting process based on entropy threshold.

        Parameters
        ----------
        thresh : float
            Entropy threshold value. Defaults to 0.
        earlystopping : bool, optional
            Whether to stop optimization early if threshold is reached. Defaults to True.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, int]
            Tuple containing arrays of Nps, entropies, and optimal Np value.
        """
        Nps = np.logspace(0, np.log10(len(self.weights)) + 1, 50, dtype=np.int64)
        entropies = np.zeros(len(Nps))
        for i in range(len(Nps)):
            self.unweight(Nps[i])
            entropies[i] = self.entropy(self.unprobs, self.reprobs)
            if entropies[i] <= thresh and earlystopping:
                loc = i
                break

        if i == len(Nps) - 1:
            try:
                loc = np.where(entropies <= thresh)[0][0]
            except:
                print("Failed minimisation procedure! Defaulting to lowest entropy.")
                loc = -1

        Nopt = Nps[loc]

        return Nps, entropies, Nopt

    def plot_entropy(self, Neff: int) -> None:
        """
        Plot the entropy as a function of the new number of replicas.
        Parameters
        ----------
        Neff : int
            Number of effective replicas
        """
        N, E, _ = self.iterate()
        fig = plt.figure()
        ax = plt.axes(xscale="log")
        ax.axvline(Neff, c="r", linestyle=":")
        ax.plot(N, E)
        ax.set_xlabel(r"Replica Number $N'_{rep}$", size=18)
        ax.set_ylabel(r"Entropy $H$", size=18)
        ax.tick_params(axis='x', direction='in', bottom=True, top=True)
        ax.tick_params(axis='y', direction='in', left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='in', bottom=True, top=True)
        ax.tick_params(axis='y', which='minor', direction='in', left=True, right=True)
        fig.savefig("Entropy.jpg", bbox_inches='tight', pad_inches=0.2)


def main(chi2: np.ndarray, N: int, store: bool = True, plot_entropy: bool = False) -> pd.DataFrame:
    """
    Perform the unweighting process and store the results.

    Parameters
    ----------
    chi2 : np.ndarray
        Array of chi-squared values.
    N : int
        Number of experimental data points that the chi2 is based on.
    store : bool, optional
        Whether to store the resulting weights in a CSV file. Defaults to True.
    plot_entropy: bool, optional
        Whether to plot and save the entropy as a function of the number of new replicas. Defaults to false.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the unweighted and reweighted values.
    """
    u = Unweight(None, (N, np.sqrt(chi2)))
    u.reweight()
    Neff = u.effective_replicas(u.reweights)
    u.unweight(Neff)

    weights = pd.DataFrame()
    weights["unweight"] = u.unweights
    weights["reweight"] = u.reweights
    weights["nrep"] = np.arange(1, len(weights) + 1)
    weights = weights.set_index("nrep", drop=True)

    if store:
        weights.to_csv("weights.csv")

    if plot_entropy:
        u.plot_entropy(Neff)

    return weights


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unweighting using chi squared values.")
    parser.add_argument("chi2_name", help="Add the filename of the chi2 dataset (.csv)")
    parser.add_argument(
        "N", help="Add the amount of experimental datapoints that the chi2 is based on"
    )
    parser.add_argument(
        "--plot_entropy", action="store_true", help="Call flag to enable entropy plotting."
    )
    args = parser.parse_args()
    chi2 = pd.read_csv(args.chi2_name).values

    main(chi2, args.N, plot_entropy=args.plot_entropy)
