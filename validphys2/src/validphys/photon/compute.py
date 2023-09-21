"""Script that calls fiatlux to add the photon PDF."""
import logging
import tempfile

import fiatlux
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.interpolate import interp1d
import yaml

from eko import beta
from eko.couplings import compute_matching_coeffs_down, compute_matching_coeffs_up
from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.n3fit_data import replica_luxseed

from . import structure_functions as sf

log = logging.getLogger(__name__)

# not the complete fiatlux runcard since some parameters are set in the code
FIATLUX_DEFAULT = {
    "apfel": False,
    "eps_rel": 1e-1,  # extra precision on any single integration.
    "mum_proton": 2.792847356,  # proton magnetic moment, from
    # http://pdglive.lbl.gov/DataBlock.action?node=S016MM which itself
    # gets it from arXiv:1203.5425 (CODATA)
    # the elastic param type, options:
    # dipole
    # A1_world_spline
    # A1_world_pol_spline
    "elastic_param": "A1_world_pol_spline",
    "elastic_electric_rescale": 1,
    "elastic_magnetic_rescale": 1,
    # the inelastic param type, options:
    "inelastic_param": "LHAPDF_Hermes_ALLM_CLAS",  # Hermes_ALLM_CLAS, LHAPDF_Hermes_ALLM_CLAS
    "rescale_r_twist4": 0,
    "rescale_r": 1,
    "allm_limits": 0,
    "rescale_non_resonance": 1,
    "rescale_resonance": 1,
    "use_mu2_as_upper_limit": False,
    "q2min_inel_override": 0.0,
    "q2max_inel_override": 1e300,
    "lhapdf_transition_q2": 9,
    # general
    "verbose": False,
}


class Photon:
    """Photon class computing the photon array with the LuxQED approach."""

    def __init__(self, theoryid, lux_params, replicas):
        theory = theoryid.get_description()
        fiatlux_runcard = FIATLUX_DEFAULT
        fiatlux_runcard["qed_running"] = bool(np.isclose(theory["Qedref"], theory["Qref"]))
        # cast explicitly from np.bool_ to bool otherwise problems in dumping it
        # TODO: for the time being, we trigger alphaem running if Qedref=Qref.
        # This is going to be changed in favor of a bool em_running
        # in the runcard
        fiatlux_runcard["mproton"] = theory["MP"]

        # precision on final integration of double integral
        if "eps_base" in lux_params:
            fiatlux_runcard["eps_base"] = lux_params["eps_base"]
            log.warning(f"Using fiatlux parameter eps_base from runcard")
        else:
            fiatlux_runcard["eps_base"] = 1e-5
            log.info(f"Using default value for fiatlux parameter eps_base")

        self.replicas = replicas

        # structure functions
        self.luxpdfset = lux_params["luxset"].load()
        self.additional_errors = lux_params["additional_errors"]
        self.luxseed = lux_params["luxseed"]

        path_to_F2 = theoryid.path / "fastkernel/FIATLUX_DIS_F2.pineappl.lz4"
        path_to_FL = theoryid.path / "fastkernel/FIATLUX_DIS_FL.pineappl.lz4"
        self.path_to_eko_photon = theoryid.path / "eko_photon.tar"
        with EKO.read(self.path_to_eko_photon) as eko:
            self.q_in = np.sqrt(eko.mu20)

        # set fiatlux
        self.lux = {}

        mb_thr = theory["kbThr"] * theory["mb"]
        mt_thr = theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100

        self.interpolator = []
        self.integral = []

        for replica in replicas:
            f2 = sf.InterpStructureFunction(path_to_F2, self.luxpdfset.members[replica])
            fl = sf.InterpStructureFunction(path_to_FL, self.luxpdfset.members[replica])
            if not np.isclose(f2.q2_max, fl.q2_max):
                log.error(
                    "FKtables for FIATLUX_DIS_F2 and FIATLUX_DIS_FL have two different q2_max"
                )

            fiatlux_runcard["q2_max"] = float(f2.q2_max)
            alpha = Alpha(theory, fiatlux_runcard["q2_max"])
            f2lo = sf.F2LO(self.luxpdfset.members[replica], theory)
            with tempfile.NamedTemporaryFile(mode="w") as tmp:
                with tmp.file as tmp_file:
                    tmp_file.write(yaml.dump(fiatlux_runcard))
                self.lux[replica] = fiatlux.FiatLux(tmp.name)
            # we have a dict but fiatlux wants a yaml file
            # TODO : once that fiatlux will allow dictionaries
            # pass directly fiatlux_runcard

            self.lux[replica].PlugAlphaQED(alpha.alpha_em, alpha.qref)
            self.lux[replica].InsertInelasticSplitQ(
                [
                    mb_thr,
                    mt_thr,
                ]
            )
            self.lux[replica].PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)

            photon_array = self.compute_photon_array(replica)
            self.interpolator.append(
                interp1d(XGRID, photon_array, fill_value="extrapolate", kind="cubic")
            )
            self.integral.append(trapezoid(photon_array, XGRID))

    def compute_photon_array(self, replica):
        r"""
        Compute the photon PDF for every point in the grid xgrid.

        Parameters
        ----------
        replica: int
            replica id

        Returns
        -------
        compute_photon_array: numpy.array
            photon PDF at the fitting scale Q0
        """
        # Compute photon PDF
        log.info(f"Computing photon")
        photon_qin = np.array(
            [self.lux[replica].EvaluatePhoton(x, self.q_in**2).total for x in XGRID]
        )
        photon_qin += self.generate_errors(replica)
        # fiatlux computes x * gamma(x)
        photon_qin /= XGRID
        # TODO : the different x points could be even computed in parallel

        # Load eko and reshape it
        with EKO.read(self.path_to_eko_photon) as eko:
            # TODO : if the eko has not the correct grid we have to reshape it
            # it has to be done inside vp-setupfit

            # construct PDFs
            pdfs_init = np.zeros((len(eko.bases.inputpids), len(XGRID)))
            for j, pid in enumerate(eko.bases.inputpids):
                if pid == 22:
                    pdfs_init[j] = photon_qin
                    ph_id = j
                else:
                    if pid not in self.luxpdfset.flavors:
                        continue
                    pdfs_init[j] = np.array(
                        [self.luxpdfset.xfxQ(x, self.q_in, replica, pid) / x for x in XGRID]
                    )

            # Apply EKO to PDFs
            for _, elem in eko.items():
                pdfs_final = np.einsum("ajbk,bk", elem.operator, pdfs_init)

        photon_Q0 = pdfs_final[ph_id]

        # we want x * gamma(x)
        return XGRID * photon_Q0

    def __call__(self, xgrid):
        """
        Compute the photon interpolating the values of self.photon_array.

        Parameters
        ----------
        xgrid : nd.array
            array of x values with shape (1,xgrid,1)

        Returns
        -------
        photon values : nd.array
            array of photon values with shape (1,xgrid,1)
        """
        return [
            self.interpolator[id](xgrid[0, :, 0])[np.newaxis, :, np.newaxis]
            for id in range(len(self.replicas))
        ]

    @property
    def error_matrix(self):
        """Generate error matrix to be used in generate_errors."""
        if not self.additional_errors:
            return None
        extra_set = self.additional_errors.load()
        qs = [self.q_in] * len(XGRID)
        res_central = np.array(extra_set.central_member.xfxQ(22, XGRID, qs))
        res = []
        for idx_member in range(101, 107 + 1):
            tmp = np.array(extra_set.members[idx_member].xfxQ(22, XGRID, qs))
            res.append(tmp - res_central)
        # first index must be x, while second one must be replica index
        return np.stack(res, axis=1)

    def generate_errors(self, replica):
        """
        Generate LUX additional errors according to the procedure
        described in sec. 2.5 of https://arxiv.org/pdf/1712.07053.pdf
        """
        if self.error_matrix is None:
            return np.zeros_like(XGRID)
        log.info(f"Generating photon additional errors")
        seed = replica_luxseed(replica, self.luxseed)
        rng = np.random.default_rng(seed=seed)
        u, s, _ = np.linalg.svd(self.error_matrix, full_matrices=False)
        errors = u @ (s * rng.normal(size=7))
        return errors


class Alpha:
    def __init__(self, theory, q2max):
        self.theory = theory
        self.alpha_s_ref = theory["alphas"]
        self.alpha_em_ref = theory["alphaqed"]
        self.qref = self.theory["Qref"]
        self.betas_qcd, self.betas_qed, self.beta_mix_qcd, self.beta_mix_qed = self.compute_betas()

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
            self.couplings_fixed_flavor = self.couplings_fixed_flavor_trn
            self.thresh, self.couplings_thresh = self.compute_couplings_at_thresholds()
        elif self.theory["ModEv"] == "EXA":
            self.couplings_fixed_flavor = self.couplings_fixed_flavor_exa
            self.thresh, self.couplings_thresh = self.compute_couplings_at_thresholds()

            xmin = XGRID[0]
            qmin = xmin * theory["MP"] / np.sqrt(1 - xmin)
            # use a lot of interpolation points since it is a long path 1e-9 -> 1e4
            self.q = np.geomspace(qmin, np.sqrt(q2max), 500, endpoint=True)

            # add threshold points in the q list since alpha is not smooth there
            self.q = np.append(self.q, [self.thresh_c, self.thresh_b, self.thresh_t])
            self.q = self.q[np.isfinite(self.q)]
            self.q.sort()

            self.alpha_vec = np.array([self.couplings_variable_flavor(q_)[1] for q_ in self.q])
            self.alpha_em = self.interpolate_alphaem
        else:
            raise ValueError(f"Evolution mode not recognized: {self.theory['ModEv']}")

    def alpha_em(self, q):
        r"""
        Compute the value of alpha_em.

        Parameters
        ----------
        q: float
            value in which alpha_em is computed

        Returns
        -------
        alpha_em: float
            electromagnetic coupling
        """
        return self.couplings_variable_flavor(q)[1]

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

    def couplings_variable_flavor(self, q):
        r"""
        Compute the value of the running couplings.

        Parameters
        ----------
        q: float
            value in which the couplings are computed

        Returns
        -------
        alpha_em: numpy.ndarray
            electromagnetic coupling
        """
        if q < self.thresh_c:
            nf = 3
        elif q < self.thresh_b:
            nf = 4
        elif q < self.thresh_t:
            nf = 5
        else:
            nf = 6
        return self.couplings_fixed_flavor(q, self.couplings_thresh[nf], self.thresh[nf], nf)

    def couplings_fixed_flavor_trn(self, q, couplings_ref, qref, nf):
        """
        Compute the running couplings for nf fixed at NLO, using truncated method.
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

        Returns
        -------
        alpha_em at NLO : float
            target value of a
        """
        alpha_ref = couplings_ref[1]
        lmu = 2 * np.log(q / qref)
        den = 1.0 + self.betas_qed[nf][0] * alpha_ref * lmu
        alpha_LO = alpha_ref / den
        alpha_NLO = alpha_LO * (1 - self.betas_qed[nf][1] * alpha_LO * np.log(den))
        return [couplings_ref[0], alpha_NLO]

    def couplings_fixed_flavor_exa(self, q, couplings_ref, qref, nf):
        """
        Compute numerically the running running for nf fixed.
        The RGEs for the two couplings are coupled by the mixed terms,
        so must be solved together starting from the same qref.

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

        Returns
        -------
        alpha_em at NLO : float
            target value of a
        """
        u = 2 * np.log(q / qref)

        # solve RGE
        res = solve_ivp(
            rge,
            (0, u),
            couplings_ref,
            args=[
                self.betas_qcd[nf],
                self.beta_mix_qcd[nf],
                self.betas_qed[nf],
                self.beta_mix_qed[nf],
            ],
            method="Radau",
            rtol=1e-6,
        )
        return [res.y[0][-1], res.y[1][-1]]

    def compute_couplings_at_thresholds(self):
        """
        Compute and store the couplings at thresholds to speed up the calling
        to alpha_em inside fiatlux:
        when q is in a certain range (e.g. thresh_c < q < thresh_b) and qref in a different one
        (e.g. thresh_b < q < thresh_t) we need to evolve from qref to thresh_b with nf=5 and then
        from thresh_b to q with nf=4. Given that the value of alpha at thresh_b is always the same
        we can avoid computing the first step storing the values of alpha in the threshold points.
        It is done for qref in a generic range (not necessarly qref=91.2).

        """
        # determine nfref
        if self.qref < self.thresh_c:
            nfref = 3
        elif self.qref < self.thresh_b:
            nfref = 4
        elif self.qref < self.thresh_t:
            nfref = 5
        else:
            nfref = 6

        thresh_list = [self.thresh_c, self.thresh_b, self.thresh_t]
        thresh_list.insert(nfref - 3, self.qref)

        thresh = {nf: thresh_list[nf - 3] for nf in range(3, self.theory["MaxNfAs"] + 1)}

        couplings_thresh = {nfref: [self.alpha_s_ref, self.alpha_em_ref]}

        logs = {
            4: np.log(self.theory["kcThr"] ** 2),
            5: np.log(self.theory["kbThr"] ** 2),
            6: np.log(self.theory["ktThr"] ** 2),
        }

        # determine the values of the couplings in the threshold points, depending on the value of qref
        for nf in range(nfref + 1, self.theory["MaxNfAs"] + 1):
            coupl_tmp = self.couplings_fixed_flavor(
                thresh[nf], couplings_thresh[nf - 1], thresh[nf - 1], nf - 1
            )
            # compute matchings
            match_up = compute_matching_coeffs_up(self.theory["HQ"], nf)
            alpha_s = coupl_tmp[0]
            alphas_new = apply_match(alpha_s, self.theory['PTO'] + 1, logs[nf], match_up)
            couplings_thresh[nf] = [alphas_new, coupl_tmp[1]]

        for nf in reversed(range(3, nfref)):
            coupl_tmp = self.couplings_fixed_flavor(
                thresh[nf], couplings_thresh[nf + 1], thresh[nf + 1], nf + 1
            )
            # compute matchings
            match_down = compute_matching_coeffs_down(self.theory["HQ"], nf - 1)
            alpha_s = coupl_tmp[0]
            alphas_new = apply_match(alpha_s, self.theory['PTO'] + 1, logs[nf + 1], match_down)
            couplings_thresh[nf] = [alphas_new, coupl_tmp[1]]

        return thresh, couplings_thresh

    def compute_betas(self):
        """Set values of betaQCD and betaQED."""
        betas_qcd = {}
        betas_qed = {}
        beta_mix_qcd = {}
        beta_mix_qed = {}
        for nf in range(3, 6 + 1):
            vec_qcd = [beta.beta_qcd_as2(nf) / (4 * np.pi)]
            vec_qed = [beta.beta_qed_aem2(nf) / (4 * np.pi)]
            for ord in range(1, self.theory['PTO'] + 1):
                vec_qcd.append(beta.b_qcd((ord + 2, 0), nf) / (4 * np.pi) ** ord)
            for ord in range(1, self.theory['QED']):
                vec_qed.append(beta.b_qed((0, ord + 2), nf) / (4 * np.pi) ** ord)
            betas_qcd[nf] = vec_qcd
            betas_qed[nf] = vec_qed
            beta_mix_qcd[nf] = beta.b_qcd((2, 1), nf) / (4 * np.pi)
            beta_mix_qed[nf] = beta.b_qed((1, 2), nf) / (4 * np.pi)
        return betas_qcd, betas_qed, beta_mix_qcd, beta_mix_qed


def apply_match(alphas, ord, L, match):
    """
    Applying matching conditions to alphas.

    Parameters
    ----------
    alphas: float
        strong coupling
    ord: int
        perturbative order
    L: float
        log(mu^2/m^2)
    match: list

    Returns
    -------
    alphas: float
        alphas after the matching

    """
    fact = 1.0
    for n in range(1, ord):
        for l_pow in range(n + 1):
            fact += (alphas / (4 * np.pi)) ** n * L**l_pow * match[n, l_pow]
    return alphas * fact


def rge(_t, alpha, beta_qcd_vec, beta_qcd_mix, beta_qed_vec, beta_qed_mix):
    """
    RGEs for the couplings. See Eqs. (5-6) of arXiv:hep-ph/9803211.
    The values of the mixed values are not used in the evolution since
    we need alpha at very low scale, below the Landau pole of alpha_s.
    This makes the alpha evolution crash. For this reason we evolve alpha
    without the mixed terms, i.e. decoupling it from alpha_s.
    rge_qcd is set to zero since we don't want an alpha_s value that is crap
    even if it is not used.
    We left the betaQCD and beta_mix part implemented in the case we find a
    solution. Anyway, it is not slowing down the code.

    """
    rge_qcd = (
        -(alpha[0] ** 2)
        * beta_qcd_vec[0]
        * (
            1
            + np.sum([alpha[0] ** (k + 1) * b for k, b in enumerate(beta_qcd_vec[1:])])
            # + alpha[1] * beta_qcd_mix
        )
    ) * 0.0
    rge_qed = (
        -(alpha[1] ** 2)
        * beta_qed_vec[0]
        * (
            1
            + np.sum([alpha[1] ** (k + 1) * b for k, b in enumerate(beta_qed_vec[1:])])
            # + alpha[0] * beta_qed_mix
        )
    )
    res = np.array([rge_qcd, rge_qed])
    return res
