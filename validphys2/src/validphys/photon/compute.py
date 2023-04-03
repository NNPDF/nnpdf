"""Script that calls fiatlux to add the photon PDF."""
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
from os import remove
import time

import fiatlux
import yaml
from eko.io import EKO
# from eko.io.manipulate import xgrid_reshape
# from eko.interpolation import XGrid

from validphys.lhapdfset import LHAPDFSet
from validphys.n3fit_data import replica_nnseed
from . import structure_functions as sf
import logging
from n3fit.io.writer import XGRID



log = logging.getLogger(__name__)

class Photon:
    """Photon class computing the photon array with the LuxQED approach."""
    def __init__(self, theoryid, fiatlux_runcard, replicas_id):
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        self.replicas_id = replicas_id
        self.q_in = 100
        self.q_in2 = self.q_in**2

        # parameters for the alphaem running
        self.alpha_em_ref = self.theory["alphaqed"]
        self.qref = self.theory["Qref"]

        self.set_betas()
        self.set_thresholds_alpha_em()

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
        for id in replicas_id:
            f2[id] = sf.InterpStructureFunction(path_to_F2, self.qcd_pdfs.members[id])
            fl[id] = sf.InterpStructureFunction(path_to_FL, self.qcd_pdfs.members[id])
            f2lo[id] = sf.F2LO(self.qcd_pdfs.members[id], self.theory)
            ff = open(f'fiatlux_runcard_{id}.yml', 'w+')
            yaml.dump(self.fiatlux_runcard, ff)
            self.lux[id] = fiatlux.FiatLux(f'fiatlux_runcard_{id}.yml')
            remove(f'fiatlux_runcard_{id}.yml')
        # we have a dict but fiatlux wants a yaml file
        # TODO : remove this dirty trick
        # we print different runcards for every replica so they do not interfere with each other
        for id in replicas_id :
            self.lux[id].PlugAlphaQED(self.alpha_em, self.qref)
            self.lux[id].InsertInelasticSplitQ([self.thresh_b, self.thresh_t if self.theory["MaxNfPdf"]==6 else 1e100])        
            self.lux[id].PlugStructureFunctions(f2[id].fxq, fl[id].fxq, f2lo[id].fxq)
        
        self.xgrid = XGRID
        self.error_matrix = self.generate_error_matrix()

        self.produce_interpolators()
    
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
        if q < self.thresh_c :
            nf = 3
        elif q < self.thresh_b :
            nf = 4
        elif q < self.thresh_t :
            nf = 5
        else :
            nf = 6
        return self.alpha_em_nlo(
            q,
            self.alpha_thresh[nf],
            self.thresh[nf],
            nf
        )
    
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
        as_NLO : float
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
        if self.theory["MaxNfAs"] <= 5 :
            self.thresh_t = np.inf
        if self.theory["MaxNfAs"] <= 4 :
            self.thresh_b = np.inf
        if self.theory["MaxNfAs"] <= 3 :
            self.thresh_c = np.inf

        thresh_list = [self.thresh_c, self.thresh_b, self.thresh_t]
        # determine nfref
        if self.qref < self.thresh_c:
            nfref = 3
        elif self.qref < self.thresh_b :
            nfref = 4
        elif self.qref < self.thresh_t :
            nfref = 5
        else :
            nfref = 6
        thresh_list.insert(nfref - 3, self.qref)

        self.thresh = {nf: thresh_list[nf - 3] for nf in range(3, self.theory["MaxNfAs"] + 1)}

        self.alpha_thresh = {nfref: self.alpha_em_ref}

        # determine the values of alpha in the threshold points, depending on the value of qref
        for nf in range(nfref + 1, self.theory["MaxNfAs"] + 1):
            self.alpha_thresh[nf] = self.alpha_em_nlo(self.thresh[nf], self.alpha_thresh[nf - 1], self.thresh[nf - 1], nf - 1)
        
        for nf in reversed(range(3, nfref)):
            self.alpha_thresh[nf] = self.alpha_em_nlo(self.thresh[nf], self.alpha_thresh[nf + 1], self.thresh[nf + 1], nf + 1)
    
    def set_betas(self):
        """Compute and store beta0 / 4pi and b1 = (beta1/beta0)/4pi as a function of nf."""
        nl = 3
        nc = 3
        eu2 = 4. / 9
        ed2 = 1. / 9
        self.beta0 = {}
        self.b1 = {}
        for nf in range(3, 6+1):
            nu = nf // 2
            nd = nf - nu
            self.beta0[nf] = ( -4.0 / 3 * (nl + nc * (nu * eu2 + nd * ed2)) ) / (4 * np.pi)
            self.b1[nf] = -4.0 * ( nl + nc * (nu * eu2**2 + nd * ed2**2) ) / self.beta0[nf] / (4 * np.pi)**2

    def compute_photon_array(self, id):
        r"""
        Compute the photon PDF for every point in the grid xgrid.

        Parameters
        ----------
        id: int
            replica id

        Returns
        -------
        compute_photon_array: numpy.array
            photon PDF at the scale 1 GeV
        """
        # Compute photon PDF
        start_time = time.perf_counter()
        photon_100GeV = np.array(
            [self.lux[id].EvaluatePhoton(x, self.q_in2).total for x in self.xgrid]
        )
        photon_100GeV += self.generate_errors(id)
        photon_100GeV /= self.xgrid
        log.info(f"Time to compute photon: {time.perf_counter() - start_time}")
        # TODO : the different x points could be even computed in parallel

        # Load eko and reshape it
        with EKO.read(self.path_to_eko_photon) as eko:
            # If we make sure that the grid of the precomputed EKO is the same of 
            # self.xgrid then we don't need to reshape
            # TODO : move the reshape inside vp-setupfit
            # xgrid_reshape(eko, targetgrid = XGrid(self.xgrid), inputgrid = XGrid(self.xgrid))

            # construct PDFs
            pdfs = np.zeros((len(eko.rotations.inputpids), len(self.xgrid)))
            for j, pid in enumerate(eko.rotations.inputpids):
                if pid == 22 :
                    pdfs[j] = photon_100GeV
                    ph_id = j
                if pid not in self.qcd_pdfs.flavors:
                    continue
                pdfs[j] = np.array(
                    [
                        self.qcd_pdfs.xfxQ(x, self.q_in, id, pid) / x
                        for x in self.xgrid
                    ]
                )
            
            # Apply EKO to PDFs
            q2 = eko.mu2grid[0]
            with eko.operator(q2) as elem:
                pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
                # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

        photon_Q0 = pdf_final[ph_id]

        # we want x * gamma(x)
        return self.xgrid * photon_Q0
    
    def produce_interpolators(self):
        """Produce the interpolation functions to be called in compute."""
        self.photons_array = [self.compute_photon_array(id) for id in self.replicas_id]
        self.interpolator = [
            interp1d(self.xgrid, photon_array, fill_value=0., kind='cubic') for photon_array in self.photons_array
        ]
    
    def compute(self, xgrid):
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
            self.interpolator[id](xgrid[0,:,0])[np.newaxis,:,np.newaxis]
            for id in range(len(self.replicas_id))
        ]
    
    def integrate(self):
        """Compute the integral of the photon on the x range."""
        return [trapezoid(self.photons_array[id], self.xgrid) for id in range(len(self.replicas_id))]
    
    def generate_error_matrix(self):
        """generate error matrix to be used for the additional errors."""
        if not self.fiatlux_runcard["additional_errors"] :
            return None
        extra_set = LHAPDFSet("LUXqed17_plus_PDF4LHC15_nnlo_100", "replicas")
        qs = [self.q_in]*len(self.xgrid)
        res_central = np.array(extra_set.central_member.xfxQ(22, self.xgrid, qs))
        res = []
        for idx_member in range(101, 107+1):
            tmp = np.array(extra_set.members[idx_member].xfxQ(22, self.xgrid, qs))
            res.append(tmp - res_central)
        # first index must be x, while second one must be replica index
        return np.stack(res, axis=1)

    def generate_errors(self, replica_id):
        """generate LUX additional errors."""
        if self.error_matrix is None :
            return np.zeros_like(self.xgrid)
        seed = replica_luxseed(replica_id, self.fiatlux_runcard["luxseed"])
        rng = np.random.default_rng(seed=seed)
        u, s, _ = np.linalg.svd(self.error_matrix, full_matrices=False)
        errors = u @ (s * rng.normal(size=7))
        return errors
    
# TODO : move it in n3fit_data.py with the others replica_seed generators?
# or use directly replica_nnseed?
def replica_luxseed(replica, luxseed):
    return replica_nnseed(replica, luxseed)
