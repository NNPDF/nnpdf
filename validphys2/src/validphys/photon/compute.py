"""Script that calls fiatlux to add the photon PDF."""
import logging
import tempfile

import fiatlux
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
import yaml

from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.n3fit_data import replica_luxseed

from . import structure_functions as sf
from .constants import ED2, EU2, NC, NL

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

        if "component" not in lux_params:
            self.component = "total"
        else:
            self.component = lux_params["component"]

        # TODO : maybe find a different name for fiatlux_dis_F2
        path_to_F2 = theoryid.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
        path_to_FL = theoryid.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"
        self.path_to_eko_photon = theoryid.path / "eko_photon.tar"
        with EKO.read(self.path_to_eko_photon) as eko:
            self.q_in = np.sqrt(eko.mu20)

        # set fiatlux
        self.lux = {}

        alpha = Alpha(theory)
        mb_thr = theory["kbThr"] * theory["mb"]
        mt_thr = theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100

        self.interpolator = {"total": [], "elastic": [], "inelastic": [], "msbar": []}
        self.integral = []

        for replica in replicas:
            f2 = sf.InterpStructureFunction(path_to_F2, self.luxpdfset.members[replica])
            fl = sf.InterpStructureFunction(path_to_FL, self.luxpdfset.members[replica])
            if not np.isclose(f2.q2_max, fl.q2_max):
                log.error(
                    "FKtables for fiatlux_dis_F2 and fiatlux_dis_FL have two different q2_max"
                )

            fiatlux_runcard["q2_max"] = float(f2.q2_max)
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

            photon_dict = self.compute_photon_array(replica)
            for key, value in photon_dict.items():
                self.interpolator[key].append(
                    interp1d(XGRID, value, fill_value="extrapolate", kind="cubic")
                )
            self.integral.append(trapezoid(photon_dict["total"], XGRID))

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
        photon_qin = {
            "total": np.zeros_like(XGRID),
            "elastic": np.zeros_like(XGRID),
            "inelastic": np.zeros_like(XGRID),
            "msbar": np.zeros_like(XGRID),
        }
        for i, x in enumerate(XGRID):
            pht = self.lux[replica].EvaluatePhoton(x, self.q_in**2)
            photon_qin["total"][i] = pht.total
            photon_qin["elastic"][i] = pht.elastic
            photon_qin["inelastic"][i] = pht.inelastic_pf
            photon_qin["msbar"][i] = pht.msbar_pf
        # photon_qin += self.generate_errors(replica)
        # TODO : the different x points could be even computed in parallel

        # Load eko and reshape it
        photon_Q0 = {}
        with EKO.read(self.path_to_eko_photon) as eko:
            # TODO : if the eko has not the correct grid we have to reshape it
            # it has to be done inside vp-setupfit

            # loop over different components
            for key, value in photon_qin.items():
                # construct PDFs
                pdfs_init = np.zeros((len(eko.bases.inputpids), len(XGRID)))
                for j, pid in enumerate(eko.bases.inputpids):
                    if pid == 22:
                        # fiatlux computes x * gamma(x)
                        pdfs_init[j] = value / XGRID
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
                # we want x * gamma(x)
                photon_Q0[key] = pdfs_final[ph_id] * XGRID

        return photon_Q0

    def __call__(self, xgrid, total):
        """
        Compute the photon interpolating the values of self.photon_array.

        Parameters
        ----------
        xgrid: nd.array
            array of x values with shape (1,xgrid,1)
        total: bool
            True for the total component, False for the others

        Returns
        -------
        photon values : nd.array
            array of photon values with shape (1,xgrid,1)
        """
        if total:
            component = "total"
        else:
            component = self.component
        return [
            self.interpolator[component][id](xgrid[0, :, 0])[np.newaxis, :, np.newaxis]
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
    def __init__(self, theory):
        self.theory = theory
        self.alpha_em_ref = theory["alphaqed"]
        self.qref = self.theory["Qedref"]
        self.beta0, self.b1 = self.set_betas()
        self.thresh, self.alpha_thresh = self.set_thresholds_alpha_em()

    def alpha_em(self, q):
        r"""
        Compute the value of alpha_em.

        Parameters
        ----------
        q: float
            value in which the coupling is computed

        Returns
        -------
        alpha_em: float
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
        return self.alpha_em_fixed_flavor(q, self.alpha_thresh[nf], self.thresh[nf], nf)

    def alpha_em_fixed_flavor(self, q, alpha_ref, qref, nf):
        """
        Compute the alpha_em running for nf fixed at NLO.

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
        lmu = 2 * np.log(q / qref)
        den = 1.0 + self.beta0[nf] * alpha_ref * lmu
        alpha_LO = alpha_ref / den
        alpha_NLO = alpha_LO * (1 - self.b1[nf] * alpha_LO * np.log(den))
        return alpha_NLO

    def set_thresholds_alpha_em(self):
        """
        Compute and store the couplings at thresholds to speed up the calling
        to alpha_em inside fiatlux:
        when q is in a certain range (e.g. thresh_c < q < thresh_b) and qref in a different one
        (e.g. thresh_b < q < thresh_t) we need to evolve from qref to thresh_b with nf=5 and then
        from thresh_b to q with nf=4. Given that the value of alpha at thresh_b is always the same
        we can avoid computing the first step storing the values of alpha in the threshold points.
        It is done for qref in a generic range (not necessarly qref=91.2).

        """
        self.thresh_c = self.theory["kcThr"] * self.theory["mc"]
        self.thresh_b = self.theory["kbThr"] * self.theory["mb"]
        self.thresh_t = self.theory["ktThr"] * self.theory["mt"]
        if self.theory["MaxNfAs"] <= 5:
            self.thresh_t = np.inf
        if self.theory["MaxNfAs"] <= 4:
            self.thresh_b = np.inf
        if self.theory["MaxNfAs"] <= 3:
            self.thresh_c = np.inf

        thresh_list = [self.thresh_c, self.thresh_b, self.thresh_t]
        # determine nfref
        if self.qref < self.thresh_c:
            nfref = 3
        elif self.qref < self.thresh_b:
            nfref = 4
        elif self.qref < self.thresh_t:
            nfref = 5
        else:
            nfref = 6
        thresh_list.insert(nfref - 3, self.qref)

        thresh = {nf: thresh_list[nf - 3] for nf in range(3, self.theory["MaxNfAs"] + 1)}

        alpha_thresh = {nfref: self.alpha_em_ref}

        # determine the values of alpha in the threshold points, depending on the value of qref
        for nf in range(nfref + 1, self.theory["MaxNfAs"] + 1):
            alpha_thresh[nf] = self.alpha_em_fixed_flavor(
                thresh[nf], alpha_thresh[nf - 1], thresh[nf - 1], nf - 1
            )

        for nf in reversed(range(3, nfref)):
            alpha_thresh[nf] = self.alpha_em_fixed_flavor(
                thresh[nf], alpha_thresh[nf + 1], thresh[nf + 1], nf + 1
            )

        return thresh, alpha_thresh

    def set_betas(self):
        """Compute and store beta0 / 4pi and b1 = (beta1/beta0)/4pi as a function of nf."""
        beta0 = {}
        b1 = {}
        for nf in range(3, 6 + 1):
            nu = nf // 2
            nd = nf - nu
            beta0[nf] = (-4.0 / 3 * (NL + NC * (nu * EU2 + nd * ED2))) / (4 * np.pi)
            b1[nf] = (
                -4.0 * (NL + NC * (nu * EU2**2 + nd * ED2**2)) / beta0[nf] / (4 * np.pi) ** 2
            )
        return beta0, b1
