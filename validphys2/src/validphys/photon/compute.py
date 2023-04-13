"""Script that calls fiatlux to add the photon PDF."""
import logging
import tempfile
import time

import fiatlux
import numpy as np
import yaml
from eko.io import EKO
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from validphys.lhapdfset import LHAPDFSet
from validphys.n3fit_data import replica_luxseed

from n3fit.io.writer import XGRID

from . import constants
from . import structure_functions as sf

log = logging.getLogger(__name__)

Q_IN = 100
EXTRA_SET = "LUXqed17_plus_PDF4LHC15_nnlo_100"


class Photon:
    """Photon class computing the photon array with the LuxQED approach."""

    def __init__(self, theoryid, fiatlux_runcard, replicas):
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        self.fiatlux_runcard["qed_running"] = "QrefQED" in self.theory
        self.replicas = replicas

        # structure functions
        self.qcd_pdfs = LHAPDFSet(fiatlux_runcard["pdf_name"], "replicas")

        # TODO : maybe find a different name for fiatlux_dis_F2
        path_to_F2 = theoryid.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
        path_to_FL = theoryid.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"
        self.path_to_eko_photon = theoryid.path / "eko_photon.tar"

        # set fiatlux
        self.lux = {}
        f2 = {}
        fl = {}
        f2lo = {}
        for replica in replicas:
            f2[replica] = sf.InterpStructureFunction(path_to_F2, self.qcd_pdfs.members[replica])
            fl[replica] = sf.InterpStructureFunction(path_to_FL, self.qcd_pdfs.members[replica])
            f2lo[replica] = sf.F2LO(self.qcd_pdfs.members[replica], self.theory)
            with tempfile.NamedTemporaryFile(mode="w") as tmp:
                with tmp.file as tmp_file:
                    tmp_file.write(yaml.dump(self.fiatlux_runcard))
                self.lux[replica] = fiatlux.FiatLux(tmp_file.name)
        # we have a dict but fiatlux wants a yaml file
        # TODO : once that fiatlux will allow dictionaries
        # pass directly self.fiatlux_runcard
        alpha = Alpha(self.theory)
        for replica in replicas:
            self.lux[replica].PlugAlphaQED(alpha.alpha_em, alpha.qref)
            self.lux[replica].InsertInelasticSplitQ(
                [
                    self.theory["kbThr"] * self.theory["mb"],
                    self.theory["ktThr"] * self.theory["mt"]
                    if self.theory["MaxNfPdf"] == 6
                    else 1e100,
                ]
            )
            self.lux[replica].PlugStructureFunctions(f2[replica].fxq, fl[replica].fxq, f2lo[replica].fxq)

        self.xgrid = XGRID
        self.error_matrix = self.generate_error_matrix()

        self.photons_array = [self.compute_photon_array(id) for id in self.replicas]
        self.interpolator = [
            interp1d(self.xgrid, photon_array, fill_value=0.0, kind="cubic")
            for photon_array in self.photons_array
        ]

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
        start_time = time.perf_counter()
        photon_qin = np.array(
            [self.lux[replica].EvaluatePhoton(x, Q_IN**2).total for x in self.xgrid]
        )
        log.info(f"Computation time: {time.perf_counter() - start_time}")
        photon_qin += self.generate_errors(replica)
        photon_qin /= self.xgrid
        # TODO : the different x points could be even computed in parallel

        # Load eko and reshape it
        with EKO.read(self.path_to_eko_photon) as eko:
            # TODO : if the eko has not the correct grid we have to reshape it
            # it has to be done inside vp-setupfit

            # construct PDFs
            pdfs_init = np.zeros((len(eko.rotations.inputpids), len(self.xgrid)))
            for j, pid in enumerate(eko.rotations.inputpids):
                if pid == 22:
                    pdfs_init[j] = photon_qin
                    ph_id = j
                if pid not in self.qcd_pdfs.flavors:
                    continue
                pdfs_init[j] = np.array(
                    [self.qcd_pdfs.xfxQ(x, Q_IN, replica, pid) / x for x in self.xgrid]
                )

            # Apply EKO to PDFs
            q2 = eko.mu2grid[0]
            with eko.operator(q2) as elem:
                pdfs_final = np.einsum("ajbk,bk", elem.operator, pdfs_init)

        photon_Q0 = pdfs_final[ph_id]

        # we want x * gamma(x)
        return self.xgrid * photon_Q0

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
            for id in range(len(self.replicas_id))
        ]

    def integrate(self):
        """Compute the integral of the photon on the x range."""
        return [
            trapezoid(self.photons_array[id], self.xgrid)
            for id in range(len(self.replicas_id))
        ]

    def generate_error_matrix(self):
        """generate error matrix to be used for the additional errors."""
        if not self.fiatlux_runcard["additional_errors"]:
            return None
        extra_set = LHAPDFSet(EXTRA_SET, "replicas")
        qs = [Q_IN] * len(self.xgrid)
        res_central = np.array(extra_set.central_member.xfxQ(22, self.xgrid, qs))
        res = []
        for idx_member in range(101, 107 + 1):
            tmp = np.array(extra_set.members[idx_member].xfxQ(22, self.xgrid, qs))
            res.append(tmp - res_central)
        # first index must be x, while second one must be replica index
        return np.stack(res, axis=1)

    def generate_errors(self, replica_id):
        """generate LUX additional errors."""
        log.info(f"Generating photon additional errors")
        if self.error_matrix is None:
            return np.zeros_like(self.xgrid)
        seed = replica_luxseed(replica_id, self.fiatlux_runcard["luxseed"])
        rng = np.random.default_rng(seed=seed)
        u, s, _ = np.linalg.svd(self.error_matrix, full_matrices=False)
        errors = u @ (s * rng.normal(size=7))
        return errors


class Alpha:
    def __init__(self, theory):
        # parameters for the alphaem running
        self.theory = theory
        self.alpha_em_ref = theory["alphaqed"]
        self.qref = self.theory.get("QrefQED", theory["Qref"])

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
        return self.alpha_em_nlo(q, self.alpha_thresh[nf], self.thresh[nf], nf)

    def alpha_em_nlo(self, q, alpha_ref, qref, nf):
        """
        Compute the alpha_em running for FFS at NLO.

        Parameters
        ----------
        q : float
            target scale
        a_ref : float
            reference value of a = alpha_em/(4*pi)
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

        thresh = {
            nf: thresh_list[nf - 3] for nf in range(3, self.theory["MaxNfAs"] + 1)
        }

        alpha_thresh = {nfref: self.alpha_em_ref}

        # determine the values of alpha in the threshold points, depending on the value of qref
        for nf in range(nfref + 1, self.theory["MaxNfAs"] + 1):
            alpha_thresh[nf] = self.alpha_em_nlo(
                thresh[nf], alpha_thresh[nf - 1], thresh[nf - 1], nf - 1
            )

        for nf in reversed(range(3, nfref)):
            alpha_thresh[nf] = self.alpha_em_nlo(
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
            beta0[nf] = (
                -4.0
                / 3
                * (
                    constants.NL
                    + constants.NC * (nu * constants.EU2 + nd * constants.ED2)
                )
            ) / (4 * np.pi)
            b1[nf] = (
                -4.0
                * (
                    constants.NL
                    + constants.NC * (nu * constants.EU2**2 + nd * constants.ED2**2)
                )
                / beta0[nf]
                / (4 * np.pi) ** 2
            )
        return beta0, b1
