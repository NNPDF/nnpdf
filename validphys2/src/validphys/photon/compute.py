"""Module that calls fiatlux to add the photon PDF."""

import logging
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
import yaml

from eko import basis_rotation
from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.n3fit_data import replica_luxseed
from validphys.loader import Loader, PhotonQEDNotFound

from . import structure_functions as sf
from .alpha import Alpha

log = logging.getLogger(__name__)
loader = Loader()

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
    """Photon class computing the photon array with the LuxQED approach.
    
    Parameters
    ----------
    theoryid : validphys.core.TheoryIDSpec
        TheoryIDSpec object describing the theory to be used as
        specified in the runcard.
    lux_params : dict
        Dictionary containing the LuxQED parameters as specified
        in the runcard.
    replica_list: list[int], optional
        List of replica ids to be computed. If None, all replicas
        will be computed based on the luxqed pdf set.
    """
    def __init__(self, theoryid, lux_params, replicas=None, save_to_disk=False, force_computation=False):
        self.theoryid = theoryid
        self.lux_params = lux_params
        self.replicas = replicas
        self.save_to_disk = save_to_disk

        fiatlux_runcard = FIATLUX_DEFAULT
        # TODO: for the time being, Qedref=Qref and so alphaem running will always trigger
        # This may be changed in the future in favor of a bool em_running in the runcard
        fiatlux_runcard["qed_running"] = True
        fiatlux_runcard["mproton"] = float(theoryid.get_description()["MP"])

        # precision on final integration of double integral
        if "eps_base" in lux_params:
            fiatlux_runcard["eps_base"] = lux_params["eps_base"]
            log.warning(f"Using fiatlux parameter eps_base from runcard")
        else:
            fiatlux_runcard["eps_base"] = 1e-5
            log.info(f"Using default value for fiatlux parameter eps_base")

        self.fiatlux_runcard = fiatlux_runcard
        # Metadata for the photon ste
        self.luxpdfset = lux_params["luxset"].load()
        self.additional_errors = lux_params["additional_errors"]
        self.luxseed = lux_params["luxseed"]
        self.luxpdfset_members = self.luxpdfset.n_members - 1 # Remove replica 0

        if force_computation:
            self.compute_photon_set()
            return
            

        try:
          self.load_photon()
        except PhotonQEDNotFound:
          log.info(f"Photon set for theory ID {self.theoryid.id} and luxset {self.luxpdfset._name} not found. Computing it now...")
          self.compute_photon_set()

    def compute_photon_set(self):
        """Compute the photon set for the desired replicas."""
        # load fiatlux
        try:
          import fiatlux
        except ModuleNotFoundError as e:
          log.error("fiatlux not found, please install fiatlux")
          raise ModuleNotFoundError("Please install fiatlux: `pip install nnpdf[qed]` or `pip install fiatlux`") from e
        
        theory = self.theoryid.get_description()

        if theory["PTO"] > 0:
            path_to_F2 = self.theoryid.path / "fastkernel/FIATLUX_DIS_F2.pineappl.lz4"
            path_to_FL = self.theoryid.path / "fastkernel/FIATLUX_DIS_FL.pineappl.lz4"

        self.path_to_eko_photon = self.theoryid.path / "eko_photon.tar"
        with EKO.read(self.path_to_eko_photon) as eko:
            self.q_in = np.sqrt(eko.mu20)

        # set fiatlux
        mb_thr = theory["kbThr"] * theory["mb"]
        mt_thr = theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100
        interpolator = []
        integral = []

        for replica in self.replicas:
            # As input replica for the photon computation we take the MOD of the luxset_members to
            # avoid failing due to limited number of replicas in the luxset
            photonreplica = (replica % self.luxpdfset_members) or self.luxpdfset_members

            f2lo = sf.F2LO(self.luxpdfset.members[photonreplica], theory)

            if theory["PTO"] > 0:
                f2 = sf.InterpStructureFunction(path_to_F2, self.luxpdfset.members[photonreplica])
                fl = sf.InterpStructureFunction(path_to_FL, self.luxpdfset.members[photonreplica])
                if not np.isclose(f2.q2_max, fl.q2_max):
                    log.error(
                        "FKtables for FIATLUX_DIS_F2 and FIATLUX_DIS_FL have two different q2_max"
                    )
                self.fiatlux_runcard["q2_max"] = float(f2.q2_max)
            else:
                f2 = f2lo
                fl = sf.FLLO()
                # using a default value for q2_max
                self.fiatlux_runcard["q2_max"] = 1e8

            alpha = Alpha(theory, self.fiatlux_runcard["q2_max"])
            with tempfile.NamedTemporaryFile(mode="w") as tmp:
                yaml.dump(self.fiatlux_runcard, tmp)
                lux = fiatlux.FiatLux(tmp.name)

            # we have a dict but fiatlux wants a yaml file
            # TODO : once that fiatlux will allow dictionaries
            # pass directly fiatlux_runcard
            lux.PlugAlphaQED(alpha.alpha_em, alpha.qref)
            lux.InsertInelasticSplitQ([mb_thr, mt_thr])
            lux.PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)

            # Evaluate photon for every point in the grid xgrid
            def evaluate_at_x(x):
              return lux.EvaluatePhoton(x, self.q_in**2).total
            
            with ThreadPoolExecutor() as executor:
                photon_qin = np.array(list(executor.map(evaluate_at_x, XGRID)))

            photon_qin += self.generate_errors(replica)

            # fiatlux computes x * gamma(x)
            photon_qin /= XGRID

            # Load eko and reshape it
            with EKO.read(self.path_to_eko_photon) as eko_photon:
                # TODO : if the eko has not the correct grid we have to reshape it
                # it has to be done inside vp-setupfit
                
                # NB: the eko should contain a single operator
                for _, elem in eko_photon.items():
                    eko_op = elem.operator

                    pdfs_init = np.zeros_like(eko_op[0, 0])
                    for j, pid in enumerate(basis_rotation.flavor_basis_pids):
                        if pid == 22:
                            pdfs_init[j] = photon_qin
                            ph_id = j
                        elif pid not in self.luxpdfset.flavors:
                            continue
                        else:
                            pdfs_init[j] = np.array(
                                [self.luxpdfset.xfxQ(x, self.q_in, photonreplica, pid) / x for x in XGRID]
                            )

                    pdfs_final = np.einsum("ajbk,bk", eko_op, pdfs_init)

            photon_Q0 = pdfs_final[ph_id]
            photon_array = XGRID * photon_Q0

            if self.save_to_disk:
                path_to_photon = loader._photons_qed_path / f"photon_theoryID_{self.theoryid.id}_fit_{self.luxpdfset._name}"
                path_to_photon.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path_to_photon / f"replica_{photonreplica}.npz",photon_array=photon_array)
                log.info(f"Saved photon replica {photonreplica} to {path_to_photon}")

            interpolator.append(interp1d(XGRID, photon_array, fill_value="extrapolate", kind="cubic"))
            integral.append(trapezoid(photon_array, XGRID))

        integral = np.stack(integral, axis=-1)


        self.integral = integral
        self.interpolator = interpolator

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
    
    def load_photon(self):
      """Load the photon resource using the Loader class."""
      path_to_photon = loader.check_photonQED(self.theoryid, self.luxpdfset._name)
      log.info(f"Loading photon QED set from {path_to_photon}")

      interpolator = []
      integral = []

      # Load the needed replicas
      for replica in self.replicas:
          # As input replica for the photon computation we take the MOD of the luxset_members to
          # avoid failing due to limited number of replicas in the luxset
          photonreplica = (replica % self.luxpdfset_members) or self.luxpdfset_members

          photon_array = np.load(path_to_photon / f"replica_{photonreplica}.npz")["photon_array"]
          interpolator.append(interp1d(XGRID, photon_array, fill_value="extrapolate", kind="cubic"))
          integral.append(trapezoid(photon_array, XGRID))
      
      self.interpolator = interpolator
      self.integral = np.stack(integral, axis=-1)
      return