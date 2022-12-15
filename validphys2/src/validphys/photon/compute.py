"""Script that calls fiatlux to add the photon PDF."""
import lhapdf
import fiatlux
import numpy as np
from eko.couplings import Couplings
from eko.compatibility import update_theory
from eko.output.legacy import load_tar
from eko.interpolation import XGrid
from eko.output.manipulate import xgrid_reshape
from eko.runner import Runner
from .structure_functions import StructureFunction

import yaml
from os import remove

class Photon:
    def __init__(self, theoryid, fiatlux_runcard, replica=0):
        # TODO : for the moment we do the 0-th replica then we change it
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        if fiatlux_runcard is not None:
            self.qref = self.theory["Qref"]
            self.q_in2 = 100**2
            theory_coupling = update_theory(self.theory)
            theory_coupling["ModEv"] = "TRN"
            theory_coupling["alphaem_running"] = True
            theory_coupling["nfref"] = 5
            theory_coupling['fact_to_ren_scale_ratio'] = 1.
            self.couplings = Couplings.from_dict(theory_coupling)

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

            self.cache = {}

    
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
        thres_charm = self.theory["Qmc"]
        thres_bottom = self.theory["Qmb"]
        thres_top = self.theory["Qmt"]
        # at LO we use ZM-VFS
        if Q < thres_charm :
            hq = 3
        elif Q < thres_bottom :
            hq = 4
        elif Q < thres_top :
            hq = 5
        else :
            hq = 6
        e2u = 4/9
        e2d = 1/9
        e2q = [e2d, e2u, e2d, e2u, e2d, e2u] # d u s c b t
        res = 0
        for i in range(1, hq+1):
            res += e2q[i-1] * (self.qcd_pdfs.xfxQ(x, Q)[i] + self.qcd_pdfs.xfxQ(x, Q)[-i])
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
        return self.couplings.a(q**2)[1] * 4 * np.pi

    def photon_fitting_scale(self, xgrids):
        r"""
        Compute the photon PDF for every point in the grid xgrid.

        Parameters
        ----------
        xgrids: numpy.array
            grid of the x points
        
        Returns
        -------
        photon_fitting_scale: numpy.array
            photon PDF at the scale 1 GeV
        """
        xgrid_list = self.exctract_grids(xgrids)
        
        grid_count = 0
        photon_list = []
        for xgrid in xgrid_list :
            photon_100GeV = np.zeros(len(xgrid))
            for i, x in enumerate(xgrid):
                for x_cached in self.cache :
                    if np.isclose(x, x_cached, atol=0, rtol=1e-3):
                        print("using cache for grid point", grid_count + i+1, "/", len(xgrids))
                        photon_100GeV[i] = self.cache[x_cached]
                        break
                else :
                    print("computing grid point", grid_count + i+1, "/", len(xgrids))
                    photon_100GeV[i] = self.lux.EvaluatePhoton(x, self.q_in2).total / x
                    self.cache[x] = photon_100GeV[i]
            grid_count += len(xgrid)

            # if the grid has less than 5 points we cannot interpolate the precomputed grid
            # since it has interpolation_polynomial_degree = 4. Therefore we compute the EKO.
            if len(xgrid) > 4:
                eko=load_tar(self.fiatlux_runcard['path_to_eko'])
                xgrid_reshape(eko, targetgrid = XGrid(xgrid), inputgrid = XGrid(xgrid))
            else :
                eko = self.compute_eko(xgrid)
            
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

            photon_fitting_scale = pdf_final[ph_id]

            # we want x * gamma(x)
            photon_list.append( xgrid * photon_fitting_scale )
       
        return np.concatenate(photon_list)
    
    def compute_eko(self, xgrid):

        theory = self.theory.copy()
        theory["nfref"] = 5
        theory["nf0"] = None
        theory["fact_to_ren_scale_ratio"] = 1.
        theory["ModSV"] = None
        theory["IB"] = 0
        theory["FNS"] = "VFNS"
        theory["QED"] = 2
        theory["ModEv"] = "EXA"
        theory["alphaem_running"] = True
        q_in = 100
        q_fin = self.theory["Q0"]
        theory["Q0"]= q_in

        operator_card = dict(
            sorted(
                dict(
                    interpolation_xgrid=xgrid.tolist(),
                    interpolation_polynomial_degree= 4 if len(xgrid) > 4 else len(xgrid) - 1,
                    interpolation_is_log=True,
                    ev_op_max_order=10,
                    ev_op_iterations=10,
                    backward_inversion="exact",
                    n_integration_cores=0,
                    debug_skip_non_singlet=False,
                    debug_skip_singlet=False,
                    Q2grid=[q_fin**2],
                    inputgrid=None,
                    targetgrid=None,
                    inputpids=None,
                    targetpids=None,
                ).items()
            )
        )

        runner = Runner(theory, operator_card)
        return runner.get_output()
    