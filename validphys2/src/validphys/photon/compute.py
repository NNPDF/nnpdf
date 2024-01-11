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
from .alpha import Alpha

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
        self.theoryid = theoryid
        self.lux_params = lux_params

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
            self.lux[replica].InsertInelasticSplitQ([mb_thr, mt_thr])
            self.lux[replica].PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)

            photon_array = self.compute_photon_array(replica)
            self.interpolator.append(
                interp1d(XGRID, photon_array, fill_value="extrapolate", kind="cubic")
            )
            self.integral.append(trapezoid(photon_array, XGRID))

        self.integral = np.stack(self.integral, axis=-1)

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
            array of photon values with shape (1, replicas, xgrid, 1)
        """
        return np.stack(
            [
                self.interpolator[id](xgrid[0, :, 0])[np.newaxis, :, np.newaxis]
                for id in range(len(self.replicas))
            ],
            axis=1,
        )

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
