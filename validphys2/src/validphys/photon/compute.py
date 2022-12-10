"""Script that calls fiatlux to add the photon PDF."""
import lhapdf
import fiatlux
import numpy as np
from eko.couplings import Couplings
from eko.compatibility import update_theory
from eko.runner import Runner
from .structure_functions import StructureFunction

import yaml

class Photon:
    def __init__(self, theoryid, fiatlux_runcard, replica=0):
        # TODO : for the moment we do the 0-th replica then we change it
        self.theory = theoryid.get_description()
        self.fiatlux_runcard = fiatlux_runcard
        if fiatlux_runcard is not None:
            # the commented objects should be passed from theory runcard,
            # however EKO complains that he doesn't find them
            self.theory["nfref"] = None
            self.theory["nf0"] = None
            self.theory["fact_to_ren_scale_ratio"] = 1.
            self.theory["ModSV"] = None
            self.theory["IC"] = 1
            self.theory["IB"] = 0
            self.theory["FNS"] = "VFNS"
            self.q_in = 100
            self.q_in2 = self.q_in ** 2
            self.q_fin = self.theory["Q0"]
            self.theory["Q0"]= self.q_in
            self.qref = self.theory["Qref"]
            self.qcd_pdfs = lhapdf.mkPDF(fiatlux_runcard["pdf_name"], replica)
            self.couplings = Couplings.from_dict(update_theory(self.theory))
            self.path_to_F2 = fiatlux_runcard["path_to_F2"]
            self.path_to_FL = fiatlux_runcard["path_to_FL"]
    
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
        mcharm = self.theory["mc"]
        mbottom = self.theory["mb"]
        mtop = self.theory["mt"]
        # at LO we use ZM-VFS
        if Q < mcharm :
            hq = 3
        elif Q < mbottom :
            hq = 4
        elif Q < mtop :
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
        xgrid: numpy.array
            grid of the x points
        pdf_name: string
            name of the QCD set
        replica: int
            number of replica
        
        Returns
        -------
        photon_fitting_scale: numpy.array
            photon PDF at the scale 1 GeV
        """
        if self.fiatlux_runcard is None :
            return np.zeros(len(xgrids))
        
        xgrid_list = self.exctract_grids(xgrids)
        
        # lux = fiatlux.FiatLux(fiatlux_runcard)
        # we have a dict but fiatlux wants a yaml file
        # TODO : remove this trick
        ff = open('fiatlux_runcard.yml', 'w+')
        yaml.dump(self.fiatlux_runcard, ff)

        lux = fiatlux.FiatLux('fiatlux_runcard.yml')
        lux.PlugAlphaQED(self.alpha_em, self.qref)

        f2 = StructureFunction(self.path_to_F2, self.qcd_pdfs)
        fl = StructureFunction(self.path_to_FL, self.qcd_pdfs)
        
        lux.PlugStructureFunctions(f2.FxQ, fl.FxQ, self.F2LO)
        
        lux.InsertInelasticSplitQ([4.18, 1e100])
        
        grid_count = 0
        photon_list = []
        for xgrid in xgrid_list :
            photon_100GeV = np.zeros(len(xgrid))
            for i, x in enumerate(xgrid):
                print("computing grid point", grid_count + i+1, "/", len(xgrids))
                photon_100GeV[i] = lux.EvaluatePhoton(x, self.q_in2).total / x
            # TODO: fiatlux returns gamma(x) or x*gamma(x) ?
            grid_count += len(xgrid)

            operator_card = dict(
                sorted(
                    dict(
                        interpolation_xgrid=xgrid.tolist(),
                        interpolation_polynomial_degree=4,
                        interpolation_is_log=True,
                        ev_op_max_order=10,
                        ev_op_iterations=10,
                        backward_inversion="expanded",
                        n_integration_cores=0,
                        debug_skip_non_singlet=False,
                        debug_skip_singlet=False,
                        Q2grid=[self.q_fin**2],
                        inputgrid=None,
                        targetgrid=None,
                        inputpids=None,
                        targetpids=None,
                    ).items()
                )
            )

            # TODO : this EKO should be precomputed and stored since it never changes
            runner = Runner(self.theory, operator_card)
            output = runner.get_output()
            
            # To be used when the new version of EKO (0.11.1) will be available on conda
            # pdfs = np.zeros((len(output.rotations.inputpids), len(xgrid)))
            # for j, pid in enumerate(output.rotations.inputpids):
            #     if pid == 22 :
            #         pdfs[j] = photon_100GeV
            #         ph_id = j
            #     if not self.qcd_pdfs.hasFlavor(pid):
            #         continue
            #     pdfs[j] = np.array(
            #         [
            #             self.qcd_pdfs.xfxQ2(pid, x, self.q_in2) / x
            #             for x in xgrid
            #         ]
            #     )
            
            pdfs = np.zeros((len(output["inputpids"]), len(output["inputgrid"])))
            for j, pid in enumerate(output["inputpids"]):
                if pid == 22 :
                    pdfs[j] = photon_100GeV
                    ph_id = j
                if not self.qcd_pdfs.hasFlavor(pid):
                    continue
                pdfs[j] = np.array(
                    [
                        self.qcd_pdfs.xfxQ2(pid, x, output["q2_ref"]) / x
                        for x in output["inputgrid"]
                    ]
                )
            
            # To be used when the new version of EKO (0.11.1) will be available on conda
            # for q2, elem in output.items():
            #     pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
            #     # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

            for q2, elem in output["Q2grid"].items():
                pdf_final = np.einsum("ajbk,bk", elem["operators"], pdfs)
                #error_final = np.einsum("ajbk,bk", elem["operator_errors"], pdfs)

            photon_fitting_scale = pdf_final[ph_id]

            # we want x * gamma(x)
            photon_list.append( xgrid * photon_fitting_scale )
        
        return np.concatenate(photon_list)
        