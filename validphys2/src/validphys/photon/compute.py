"""Module that calls fiatlux to add the photon PDF."""

import logging
import tempfile
from joblib import Parallel, delayed
from functools import lru_cache

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
import yaml

from eko import basis_rotation
from eko.io import EKO
from n3fit.io.writer import XGRID
from validphys.n3fit_data import replica_luxseed
from validphys.loader import Loader, PhotonQEDNotFound
from validphys.core import FKTableSpec

from . import structure_functions as sf
from .alpha import Alpha

log = logging.getLogger(__name__)

# not the complete fiatlux runcard since some parameters are set in the code
FIATLUX_DEFAULT = {
    "apfel": False,
    "eps_base": 1e-5,
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
    def __init__(self, theoryid, lux_params, replicas: list[int], save_to_disk=False, force_computation=False, q_in=100.0):
        self.theoryid = theoryid
        self.lux_params = lux_params
        self.replicas = replicas
        self.save_to_disk = save_to_disk
        self.q_in = q_in
        self.luxpdfset = lux_params["luxset"].load()
        self.luxpdfset_members = self.luxpdfset.n_members - 1 #

        # Compute or load photon_qin
        if force_computation:
          photon_in = self._compute_photon_set()
        else:
          try:
            photon_in = self._load_photon()
          except PhotonQEDNotFound:
            photon_in = self._compute_photon_set()
        self.photon_qin = photon_in

    def _load_photon(self):
      """Load the photon resource using the Loader class."""
      loader = Loader()
      path_to_photon = loader.check_photonQED(self.theoryid.id, self.luxpdfset._name)
      log.info(f"Loading photon QED set from {path_to_photon}")

      # Load the needed replicas
      photon_qin_array = []
      for replica in self.replicas:
          # As input replica for the photon computation we take the MOD of the luxset_members to
          # avoid failing due to limited number of replicas in the luxset
          photonreplica = (replica % self.luxpdfset_members) or self.luxpdfset_members
          photon_qin = np.load(path_to_photon / f"replica_{photonreplica}.npz")["photon_qin"]
          photon_qin_array.append(photon_qin)

      photon_qin_array = np.stack(photon_qin_array, axis=0)
      return photon_qin_array
    
    def _compute_photon_set(self):
       try:
          import fiatlux
       except ModuleNotFoundError as e:              
          log.error("fiatlux not found, please install fiatlux")
          raise ModuleNotFoundError("Please install fiatlux: `pip install nnpdf[qed]` or `pip install fiatlux`") from e
       
       #########
       replicas = self.replicas
       luxset = self.lux_params['luxset'].load()
       fiatlux_runcard = self._setup_fiatlux_runcard()
       #########

       # Set fiatlux parameters
       theory = self.theoryid.get_description()
       mb_thr = theory["kbThr"] * theory["mb"]
       mt_thr = theory["ktThr"] * theory["mt"] if theory["MaxNfPdf"] == 6 else 1e100
       f2lo_func = lambda member: sf.F2LO(member, theory)
       if theory["PTO"] > 0:
            path_to_F2 = self.theoryid.path / "fastkernel/FIATLUX_DIS_F2.pineappl.lz4"
            path_to_FL = self.theoryid.path / "fastkernel/FIATLUX_DIS_FL.pineappl.lz4"
            f2_func = lambda member: sf.InterpStructureFunction(path_to_F2, member)
            fl_func = lambda member: sf.InterpStructureFunction(path_to_FL, member)
            
            # Use central replica to check q2_max
            f2 = f2_func(luxset.central_member)
            fl = fl_func(luxset.central_member)
            if not np.isclose(f2.q2_max, fl.q2_max):
                log.error(
                    "FKtables for FIATLUX_DIS_F2 and FIATLUX_DIS_FL have two different q2_max"
                )
                raise ValueError("FKtables for FIATLUX_DIS_F2 and FIATLUX_DIS_FL have two different q2_max")
            fiatlux_runcard["q2_max"] = float(f2.q2_max)
       else:
            f2_func = f2lo_func
            fl_func = lambda _: sf.FLLO()

       alpha = Alpha(theory, fiatlux_runcard["q2_max"])

       if self.save_to_disk:
          loader = Loader()
          path_to_photon = loader._photons_qed_path / f"photon_theoryID_{self.theoryid.id}_fit_{luxset._name}"
          path_to_photon.mkdir(parents=True, exist_ok=True)

       photon_qin_array = []
       for replica in replicas:
          # Avoid failing due to limited number of replicas in the luxset
          photonreplica = (replica % self.luxpdfset_members) or self.luxpdfset_members
          f2lo = f2lo_func(luxset.members[photonreplica])
          f2 = f2_func(luxset.members[photonreplica])
          fl = fl_func(luxset.members[photonreplica])

          with tempfile.NamedTemporaryFile(mode="w") as tmp:
              yaml.dump(fiatlux_runcard, tmp)
              lux = fiatlux.FiatLux(tmp.name)

          # we have a dict but fiatlux wants a yaml file
          # TODO : once that fiatlux will allow dictionaries
          # pass directly fiatlux_runcard
          lux.PlugAlphaQED(alpha.alpha_em, alpha.qref)
          lux.InsertInelasticSplitQ([mb_thr, mt_thr])
          lux.PlugStructureFunctions(f2.fxq, fl.fxq, f2lo.fxq)

          photon_qin = np.array(
            [lux.EvaluatePhoton(x, self.q_in**2).total for x in XGRID]
          )
          photon_qin += self._generate_errors(replica)

          # fiatlux computes x * gamma(x)
          photon_qin /= XGRID

          if self.save_to_disk:
              np.savez_compressed(path_to_photon / f"replica_{photonreplica}.npz", photon_qin=photon_qin)
              log.info(f"Saved photon replica {photonreplica} to {path_to_photon}")

          photon_qin_array.append(photon_qin)

       photon_qin_array = np.stack(photon_qin_array, axis=0)
       return photon_qin_array
       
    def _setup_fiatlux_runcard(self):
      fiatlux_runcard = FIATLUX_DEFAULT

      # TODO: for the time being, Qedref=Qref and so alphaem running will always trigger
      # This may be changed in the future in favor of a bool em_running in the runcard
      fiatlux_runcard["qed_running"] = True
      fiatlux_runcard["mproton"] = float(self.theoryid.get_description()["MP"])

      # precision on final integration of double integral
      if "eps_base" in self.lux_params:
          fiatlux_runcard["eps_base"] = self.lux_params["eps_base"]
          log.warning(f"Using fiatlux parameter eps_base from runcard")
      else:
          log.info(f"Using default value for fiatlux parameter eps_base = {fiatlux_runcard['eps_base']}")
      return fiatlux_runcard
    
    @property
    def error_matrix(self):
        """Generate error matrix to be used in generate_errors."""
        if "additional_errors" not in self.lux_params:
            return None
        extra_set = self.lux_params["additional_errors"].load()
        qs = [self.q_in] * len(XGRID)
        res_central = np.array(extra_set.central_member.xfxQ(22, XGRID, qs))
        res = []
        for idx_member in range(101, 107 + 1):
            tmp = np.array(extra_set.members[idx_member].xfxQ(22, XGRID, qs))
            res.append(tmp - res_central)
        # first index must be x, while second one must be replica index
        return np.stack(res, axis=1)
    
    def _generate_errors(self, replica):
        """
        Generate LUX additional errors according to the procedure
        described in sec. 2.5 of https://arxiv.org/pdf/1712.07053.pdf
        """
        if self.error_matrix is None:
            return np.zeros_like(XGRID)
        log.info(f"Generating photon additional errors")
        luxseed = self.lux_params.get("luxseed", None)
        seed = replica_luxseed(replica, luxseed)
        rng = np.random.default_rng(seed=seed)
        u, s, _ = np.linalg.svd(self.error_matrix, full_matrices=False)
        errors = u @ (s * rng.normal(size=7))
        return errors

    @lru_cache(maxsize=None)
    def _evolve(self):
        """Perform the EKOs to evolve the photon from q_in to Q0."""
        log.info(f"Evolving photon from q_in={self.q_in} GeV to Q0 using EKO...")
        photon_qin_array = self.photon_qin
        interpolator = []
        integral = []

        path_to_eko_photon = self.theoryid.path / "eko_photon.tar"
        with EKO.read(path_to_eko_photon) as eko_photon:
            # Check that qin mathces with the one in the EKO
            if not np.isclose(self.q_in, np.sqrt(eko_photon.mu20)):
                log.error(f"Photon q_in {self.q_in} does not match the one in the EKO {np.sqrt(eko_photon.mu20)}")
                raise ValueError("Photon q_in does not match the one in the EKO")

            # TODO : if the eko has not the correct grid we have to reshape it
            # it has to be done inside vp-setupfit
            for replica in self.replicas:
              photonreplica = (replica % self.luxpdfset_members) or self.luxpdfset_members
            
              # NB: the eko should contain a single operator
              for _, elem in eko_photon.items():
                  eko_op = elem.operator

                  pdfs_init = np.zeros_like(eko_op[0, 0])
                  for j, pid in enumerate(basis_rotation.flavor_basis_pids):
                      if pid == 22:
                          pdfs_init[j] = photon_qin_array[replica-1]
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
              interpolator.append(interp1d(XGRID, photon_array, fill_value="extrapolate", kind="cubic"))
              integral.append(trapezoid(photon_array, XGRID))

        integral = np.stack(integral, axis=-1)
        return integral, interpolator
    
    @property
    def integral(self):
        """Return the integral values."""
        integral, _ = self._evolve()
        return integral

    @property
    def interpolator(self):
        """Return the interpolator functions."""
        _, interpolator = self._evolve()
        return interpolator
    
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
    
def compute_photon_to_disk(theoryid, fiatlux, maxcores):
    """Function to compute the photon PDF set.""" 
    luxset = fiatlux['luxset'].load()
    force_compute = fiatlux.get('compute_in_setupfit', False)
    replicas = list(range(1, luxset.n_members))

    # Return None and avoid pickling issues with the photon class.
    def wrapper_fn(replica, theoryid, fiatlux, force_compute):
        _ = Photon(theoryid, fiatlux, replicas=[replica], save_to_disk=force_compute, force_computation=force_compute)
        return None
    
    log.info(f"Starting computation of the photon using {maxcores} effective cores...")
    _ = Parallel(n_jobs=maxcores)(delayed(wrapper_fn)(replica=replica, theoryid=theoryid, fiatlux=fiatlux, force_compute=force_compute) for replica in replicas)