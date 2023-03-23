"""Script that calls fiatlux to add the photon PDF."""
import lhapdf

import numpy as np

from . import structure_functions as sf
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid

from eko.io import EKO
from eko.io.manipulate import xgrid_reshape
from eko.interpolation import XGrid

from n3fit.io.writer import XGRID

import fiatlux

import yaml
from os import remove
import time

class Photon:
    def __init__(self, theoryid, fiatlux_runcard, replicas_id):
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        self.replicas_id = replicas_id
        self.q_in2 = 100**2

        # parameters for the alphaem running
        self.alpha_em_ref = self.theory["alphaqed"]
        self.qref = self.theory["Qref"]
        # TODO : maybe they shoud be kDIS instead of k, but usually they are the same
        self.thresh_c = self.theory["kcThr"] * self.theory["mc"]
        self.thresh_b = self.theory["kbThr"] * self.theory["mb"]
        self.thresh_t = self.theory["ktThr"] * self.theory["mt"]
        if self.theory["MaxNfAs"] <= 5 :
            self.thresh_t = np.inf
        if self.theory["MaxNfAs"] <= 4 :
            self.thresh_b = np.inf
        if self.theory["MaxNfAs"] <= 3 :
            self.thresh_c = np.inf
        self.set_betas()
        self.set_thresholds_alpha_em()

        # structure functions
        self.qcd_pdfs = [lhapdf.mkPDF(fiatlux_runcard["pdf_name"], id) for id in replicas_id]

        # TODO : maybe find a different name for fiatlux_dis_F2
        path_to_F2 = theoryid.path / "fastkernel/fiatlux_dis_F2.pineappl.lz4"
        path_to_FL = theoryid.path / "fastkernel/fiatlux_dis_FL.pineappl.lz4"
        f2 = [sf.StructureFunction(path_to_F2, pdfs) for pdfs in self.qcd_pdfs]
        fl = [sf.StructureFunction(path_to_FL, pdfs) for pdfs in self.qcd_pdfs]
        f2lo = [sf.F2LO(pdfs, self.theory) for pdfs in self.qcd_pdfs]

        self.path_to_eko_photon = theoryid.path / "eko_photon.tar"

        # set fiatlux
        ff = open('fiatlux_runcard.yml', 'w+')
        yaml.dump(self.fiatlux_runcard, ff)
        self.lux = [fiatlux.FiatLux('fiatlux_runcard.yml') for i in range(len(replicas_id))]
        remove('fiatlux_runcard.yml')
        # we have a dict but fiatlux wants a yaml file
        # TODO : remove this dirty trick
        for i in range(len(replicas_id)):
            self.lux[i].PlugAlphaQED(self.alpha_em, self.qref)
            self.lux[i].InsertInelasticSplitQ([self.thresh_b, self.thresh_t if self.theory["MaxNfPdf"]==6 else 1e100])        
            self.lux[i].PlugStructureFunctions(f2[i].FxQ, fl[i].FxQ, f2lo[i].FxQ)
        
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
        """Compute and store the couplings at thresholds"""
        thresh_list = [self.thresh_c, self.thresh_b, self.thresh_t]
        if self.qref < self.thresh_c :
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
        xgrids: numpy.array
            grid of the x points
        
        Returns
        -------
        compute_photon_array: numpy.array
            photon PDF at the scale 1 GeV
        """
        # Compute photon PDF
        start_time = time.perf_counter()
        photon_100GeV = np.array([self.lux[id].EvaluatePhoton(x, self.q_in2).total for x in self.xgrid])
        photon_100GeV += self.generate_errors()
        photon_100GeV /= self.xgrid
        print("Time to compute photon:", time.perf_counter() - start_time)
        # TODO : the different x points could be even computed in parallel

        # Load eko and reshape it
        with EKO.edit(self.path_to_eko_photon) as eko:
            # If we make sure that the grid of the precomputed EKO is the same of 
            # self.xgrid then we don't need to reshape
            xgrid_reshape(eko, targetgrid = XGrid(self.xgrid), inputgrid = XGrid(self.xgrid))

            # construct PDFs
            pdfs = np.zeros((len(eko.rotations.inputpids), len(self.xgrid)))
            for j, pid in enumerate(eko.rotations.inputpids):
                if pid == 22 :
                    pdfs[j] = photon_100GeV
                    ph_id = j
                if not self.qcd_pdfs[id].hasFlavor(pid):
                    continue
                pdfs[j] = np.array(
                    [
                        self.qcd_pdfs[id].xfxQ2(pid, x, self.q_in2) / x
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
        self.photons_array = [self.compute_photon_array(i) for i in range(len(self.replicas_id))]
        self.interpolator = [interp1d(self.xgrid, photon_array, fill_value=0., kind='cubic') for photon_array in self.photons_array]
    
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
        extra_set = lhapdf.mkPDFs("LUXqed17_plus_PDF4LHC15_nnlo_100")
        return np.array(
            [
                [(extra_set[i].xfxQ2(22, x, self.q_in2) - extra_set[0].xfxQ2(22, x, self.q_in2)) for i in range(101, 107+1)]
                for x in self.xgrid
            ]
        ) # first index must be x, while second one must be replica index

    def generate_errors(self):
        """generate LUX additional errors."""
        if self.error_matrix is None :
            return np.zeros_like(self.xgrid)
        u, s, _ = np.linalg.svd(self.error_matrix, full_matrices=False)
        errors = u @ (s * np.random.normal(size=7))
        return errors
