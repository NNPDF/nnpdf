"""Script that calls fiatlux to add the photon PDF."""
import lhapdf
import fiatlux
import numpy as np
from eko.output.legacy import load_tar
from eko.interpolation import XGrid
from eko.output.manipulate import xgrid_reshape
from eko.interpolation import make_grid
from .structure_functions import StructureFunction
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid

import yaml
from os import remove

class Photon:
    def __init__(self, theoryid, fiatlux_runcard, replica=0):
        # TODO : for the moment we do the 0-th replica then we change it
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        if fiatlux_runcard is not None:
            self.q_in2 = 100**2
            self.alpha_em_ref = self.theory["alphaqed"]
            self.qref = self.theory["Qref"]
            self.eu2 = 4. / 9
            self.ed2 = 1. / 9
            self.e2q = [self.ed2, self.eu2, self.ed2, self.eu2, self.ed2, self.eu2] # d u s c b t
            self.set_thresholds_a_em()

            self.qcd_pdfs = lhapdf.mkPDF(fiatlux_runcard["pdf_name"], replica)
            path_to_F2 = fiatlux_runcard["path_to_F2"]
            path_to_FL = fiatlux_runcard["path_to_FL"]
            f2 = StructureFunction(path_to_F2, self.qcd_pdfs)
            fl = StructureFunction(path_to_FL, self.qcd_pdfs)

            # lux = fiatlux.FiatLux(fiatlux_runcard)
            # we have a dict but fiatlux wants a yaml file
            # TODO : remove this dirty trick
            ff = open('fiatlux_runcard.yml', 'w+')
            yaml.dump(self.fiatlux_runcard, ff)
            self.lux = fiatlux.FiatLux('fiatlux_runcard.yml')
            remove('fiatlux_runcard.yml')
            self.lux.PlugAlphaQED(self.alpha_em, self.qref)            
            self.lux.PlugStructureFunctions(f2.FxQ, fl.FxQ, self.F2LO)
            self.lux.InsertInelasticSplitQ([4.18, 1e100])

            self.produce_interpolator()

    
    def exctract_grids(self, xgrids):
        r"""
        Extract the subgrids inside xgrids.

        xgrids is the concatenation of different grids, i.e.
        xgrid = np.array([xmin1, ..., xmax1, xmin2, ...,xmax2, xmin3, ...]).
        The different grids are extracted and stored in a list:
        xgrid_list = [np.array([xgrid1]), np.array([xgrid2]), ...]

        Parameters
        ----------
        xgrids : nd.array
            concatenation of the subgrids
        
        Returns
        -------
        xgrid_list : list
            list containing the different grids
        """
        xgrid_list = []
        imin = 0
        for i in range(1, len(xgrids)):
            if xgrids[i-1] > xgrids[i] :
                xgrid_list.append(xgrids[imin:i])
                imin = i
        xgrid_list.append(xgrids[imin:])
        return xgrid_list
    
    def F2LO(self, x, Q):
        r"""
        Compute the LO DIS structure function F2.

        Parameters
        ----------
        x : float
            Bjorken's variable
        Q : float
            DIS hard scale
        
        Returns
        -------
        F2_LO : float
            Structure function F2 at LO
        """
        # at LO we use ZM-VFS
        if Q < self.theory["Qmc"] :
            nf = 3
        elif Q < self.theory["Qmb"] :
            nf = 4
        elif Q < self.theory["Qmt"] :
            nf = 5
        else :
            nf = 6
        res = 0
        for i in range(1, nf+1):
            res += self.e2q[i-1] * (self.qcd_pdfs.xfxQ(x, Q)[i] + self.qcd_pdfs.xfxQ(x, Q)[-i])
        return res

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
        # return self.couplings.a(q**2)[1] * 4 * np.pi
        if q < self.theory["Qmc"] :
            nf = 3
        elif q < self.theory["Qmb"] :
            nf = 4
        elif q < self.theory["Qmt"] :
            nf = 5
        else :
            nf = 6
        return self.a_em_nlo(
            q,
            self.a_thresh[nf],
            self.thresh[nf],
            nf
        ) * (4 * np.pi)
    
    def a_em_nlo(self, q, a_ref, qref, nf):
        nl = 3
        nc = 3
        nu = nf // 2
        nd = nf - nu
        beta0 = ( -4.0 / 3 * (nl + nc * (nu * self.eu2 + nd * self.ed2)) )
        beta1 = -4.0 * ( nl + nc * (nu * self.eu2**2 + nd * self.ed2**2) )
        lmu = np.log(q / qref)
        den = 1.0 + beta0 * a_ref * lmu
        a_LO = a_ref / den
        as_NLO = a_LO * (1 - beta1 / beta0 * a_LO * np.log(den))
        return as_NLO
    
    def set_thresholds_a_em(self):
        a_ref = self.alpha_em_ref / (4 * np.pi)
        self.a_em_mt = self.a_em_nlo(self.theory["Qmt"], a_ref, self.qref, 5)
        self.a_em_mb = self.a_em_nlo(self.theory["Qmb"], a_ref, self.qref, 5)
        self.a_em_mc = self.a_em_nlo(self.theory["Qmc"], self.a_em_mb, self.theory["Qmb"], 4)

        self.thresh = {3: self.theory["Qmc"], 4: self.theory["Qmb"], 5: self.qref, 6:self.theory["Qmt"]}
        self.a_thresh = {3: self.a_em_mc, 4:self.a_em_mb, 5:self.alpha_em_ref/(4*np.pi), 6:self.a_em_mt}

    def compute_photon_array(self, xgrids):
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
        xgrid_list = self.exctract_grids(xgrids)
        
        photon_list = []
        for xgrid in xgrid_list :
            photon_100GeV = np.zeros(len(xgrid))
            for i, x in enumerate(xgrid):
                print("computing grid point", i+1, "/", len(xgrids))
                photon_100GeV[i] = self.lux.EvaluatePhoton(x, self.q_in2).total / x

            eko=load_tar(self.fiatlux_runcard['path_to_eko'])
            xgrid_reshape(eko, targetgrid = XGrid(xgrid), inputgrid = XGrid(xgrid))
            
            pdfs = np.zeros((len(eko.rotations.inputpids), len(xgrid)))
            for j, pid in enumerate(eko.rotations.inputpids):
                if pid == 22 :
                    pdfs[j] = photon_100GeV
                    ph_id = j
                if not self.qcd_pdfs.hasFlavor(pid):
                    continue
                pdfs[j] = np.array(
                    [
                        self.qcd_pdfs.xfxQ2(pid, x, self.q_in2) / x
                        for x in xgrid
                    ]
                )
            
            for q2, elem in eko.items():
                pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
                # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

            photon_Q0 = pdf_final[ph_id]
            # we want x * gamma(x)
            photon_list.append( xgrid * photon_Q0 )
       
        return np.concatenate(photon_list)
    
    def produce_interpolator(self):
        self.xgrid = make_grid(98, 99, x_min=1.e-9) # TODO : use the output grid so the EKO will not be reshaped
        self.photon_array = self.compute_photon_array(self.xgrid)
        self.interpolator = interp1d(self.xgrid, self.photon_array, fill_value=0.)
    
    def compute(self, xgrid):
        return self.interpolator(xgrid[0,:,0])[np.newaxis,:,np.newaxis]
    
    def integrate(self):
        return trapezoid(self.photon_array, self.xgrid)