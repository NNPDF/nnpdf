"""Module that implements the alpha running"""
import numpy as np
from scipy.integrate import solve_ivp

from eko import beta
from n3fit.io.writer import XGRID

from .constants import ME, MMU, MQL, MTAU


class Alpha:
    """Class that solves the RGE for alpha_qed"""
    def __init__(self, theory, q2max):
        self.theory = theory
        self.alpha_em_ref = theory["alphaqed"]
        self.qref = theory["Qref"]
        self.qed_order = theory['QED']

        # compute and store thresholds
        self.thresh_c = self.theory["kcThr"] * self.theory["mc"]
        self.thresh_b = self.theory["kbThr"] * self.theory["mb"]
        self.thresh_t = self.theory["ktThr"] * self.theory["mt"]
        if self.theory["MaxNfAs"] <= 5:
            self.thresh_t = np.inf
        if self.theory["MaxNfAs"] <= 4:
            self.thresh_b = np.inf
        if self.theory["MaxNfAs"] <= 3:
            self.thresh_c = np.inf

        if self.theory["ModEv"] == "TRN":
            self.alphaem_fixed_flavor = self.alphaem_fixed_flavor_trn
            self.thresh, self.alphaem_thresh = self.compute_alphaem_at_thresholds()
        elif self.theory["ModEv"] == "EXA":
            self.alphaem_fixed_flavor = self.alphaem_fixed_flavor_exa
            self.thresh, self.alphaem_thresh = self.compute_alphaem_at_thresholds()

            xmin = XGRID[0]
            qmin = xmin * theory["MP"] / np.sqrt(1 - xmin)
            # use a lot of interpolation points since it is a long path 1e-9 -> 1e4
            self.q = np.geomspace(qmin, np.sqrt(q2max), 500, endpoint=True)

            # add threshold points in the q list since alpha is not smooth there
            self.q = np.append(
                self.q, [ME, MMU, MQL, self.thresh_c, MTAU, self.thresh_b, self.qref, self.thresh_t]
            )
            self.q = self.q[np.isfinite(self.q)]
            self.q.sort()

            self.alpha_vec = np.array([self.alpha_em(q_) for q_ in self.q])
            self.alpha_em = self.interpolate_alphaem
        else:
            raise ValueError(f"Evolution mode not recognized: {self.theory['ModEv']}")

    def interpolate_alphaem(self, q):
        r"""
        Interpolate precomputed values of alpha_em.

        Parameters
        ----------
        q: float
            value in which alpha_em is computed

        Returns
        -------
        alpha_em: float
            electromagnetic coupling
        """
        return np.interp(q, self.q, self.alpha_vec)

    def alpha_em(self, q):
        r"""
        Compute the value of the running alphaem.

        Parameters
        ----------
        q: float
            value in which alphaem is computed

        Returns
        -------
        alpha_em: numpy.ndarray
            electromagnetic coupling
        """
        nf, nl = self.find_region(q)
        return self.alphaem_fixed_flavor(
            q, self.alphaem_thresh[(nf, nl)], self.thresh[(nf, nl)], nf, nl
        )

    def alphaem_fixed_flavor_trn(self, q, alphaem_ref, qref, nf, nl):
        """
        Compute the running alphaem for nf fixed at, at least, NLO, using truncated method.
        In this case the RGE for alpha_em is solved decoupling it from the RGE for alpha_s
        (so the mixed terms are removed). alpha_s will just be unused.

        Parameters
        ----------
        q : float
            target scale
        alph_aem_ref : float
            reference value of alpha_em
        qref: float
            reference scale
        nf: int
            number of flavors
        nl: int
            number of leptons

        Returns
        -------
        alpha_em at NLO : float
            target value of a
        """
        if (nf, nl) == (0, 0):
            return alphaem_ref
        lmu = 2 * np.log(q / qref)
        den = 1.0 + self.betas_qed[(nf, nl)][0] * alphaem_ref * lmu
        alpha_LO = alphaem_ref / den
        alpha_NLO = alpha_LO * (
            1 - self.betas_qed[(nf, nl)][1] / self.betas_qed[(nf, nl)][0] * alpha_LO * np.log(den)
        )
        return alpha_NLO

    def alphaem_fixed_flavor_exa(self, q, alphaem_ref, qref, nf, nl):
        """
        Compute numerically the running alphaem for nf fixed.

        Parameters
        ----------
        q : float
            target scale
        alph_aem_ref : float
            reference value of alpha_em
        qref: float
            reference scale
        nf: int
            number of flavors
        nl: int
            number of leptons

        Returns
        -------
        alpha_em: float
            target value of a
        """
        # at LO in QED the TRN solution is exact
        if self.qed_order == 1:
            return self.alphaem_fixed_flavor_trn(q, alphaem_ref, qref, nf, nl)
        
        u = 2 * np.log(q / qref)

        # solve RGE
        res = solve_ivp(
            _rge, (0, u), (alphaem_ref,), args=[self.betas_qed[(nf, nl)]], method="Radau", rtol=1e-6
        )
        # taking fist (and only) element of y since it is a 1-D differential equation.
        #  then the last element of the array which corresponds to alpha(q)
        return res.y[0][-1]

    def compute_alphaem_at_thresholds(self):
        """
        Compute and store alphaem at thresholds to speed up the calling
        to alpha_em inside fiatlux:
        when q is in a certain range (e.g. thresh_c < q < thresh_b) and qref in a different one
        (e.g. thresh_b < q < thresh_t) we need to evolve from qref to thresh_b with nf=5 and then
        from thresh_b to q with nf=4. Given that the value of alpha at thresh_b is always the same
        we can avoid computing the first step storing the values of alpha in the threshold points.
        It is done for qref in a generic range (not necessarly qref=91.2).

        """
        # determine nfref
        nfref, nlref = self.find_region(self.qref)

        thresh_list = [ME, MMU, MQL, self.thresh_c, MTAU, self.thresh_b, self.thresh_t]
        thresh_list.append(self.qref)
        thresh_list.sort()

        thresh = {}

        for mq in thresh_list:
            eps = -1e-6 if mq < self.qref else 1e-6
            if np.isfinite(mq):
                nf_, nl_ = self.find_region(mq + eps)
                thresh[(nf_, nl_)] = mq

        regions = list(thresh.keys())

        self.betas_qed = self.compute_betas(regions)

        start = regions.index((nfref, nlref))

        alphaem_thresh = {(nfref, nlref): self.alpha_em_ref}

        for i in range(start + 1, len(regions)):
            alphaem_thresh[regions[i]] = self.alphaem_fixed_flavor(
                thresh[regions[i]],
                alphaem_thresh[regions[i - 1]],
                thresh[regions[i - 1]],
                *regions[i - 1],
            )

        for i in reversed(range(0, start)):
            alphaem_thresh[regions[i]] = self.alphaem_fixed_flavor(
                thresh[regions[i]],
                alphaem_thresh[regions[i + 1]],
                thresh[regions[i + 1]],
                *regions[i + 1],
            )

        return thresh, alphaem_thresh

    def compute_betas(self, regions):
        """Set values of betaQCD and betaQED."""
        betas_qed = {}
        for nf, nl in regions:
            betas_qed[(nf, nl)] = [
                beta.beta_qed_aem2(nf, nl) / (4 * np.pi),
                beta.beta_qed_aem3(nf, nl) / (4 * np.pi) ** 2 if self.theory['QED'] == 2 else 0.,
            ]
        return betas_qed

    def find_region(self, q):
        if q < ME:
            nf = 0
            nl = 0
        elif q < MMU:
            nf = 0
            nl = 1
        elif q < MQL:
            nf = 0
            nl = 2
        elif q < self.thresh_c:
            nf = 3
            nl = 2
        elif q < MTAU:
            nf = 4
            nl = 2
        elif q < self.thresh_b:
            nf = 4
            nl = 3
        elif q < self.thresh_t:
            nf = 5
            nl = 3
        else:
            nf = 6
            nl = 3
        return nf, nl


def _rge(_t, alpha_t, beta_qed_vec):
    """RGEs for the running of alphaem"""
    rge_qed = -(alpha_t**2) * (
        beta_qed_vec[0] + np.sum([alpha_t ** (k + 1) * beta for k, beta in enumerate(beta_qed_vec[1:])])
    )
    return rge_qed
