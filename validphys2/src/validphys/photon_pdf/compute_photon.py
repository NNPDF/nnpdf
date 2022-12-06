"""Script that calls fiatlux to add the photon PDF."""

# from validphys import lhapdfset
from validphys.api import API
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
        self.theory = API.theoryid(theoryid = theoryid).get_description().copy()
        self.fiatlux_runcard = fiatlux_runcard
        # theory["nfref"] = None
        # theory["nf0"] = None
        # theory["fact_to_ren_scale_ratio"] = 1.
        self.theory["ModSV"] = None
        # theory["IC"]=0
        # theory["IB"]=0
        self.theory["FNS"] = "VFNS"
        self.q_in = 100
        self.q_in2 = self.q_in ** 2
        self.q_fin = self.theory["Q0"]
        self.theory["Q0"]= self.q_in
        self.qref = self.theory["Qref"]
        # xir = theory["XIR"]
        # xif = theory["XIF"]
        self.qcd_pdfs = lhapdf.mkPDF(fiatlux_runcard["pdf_name"], replica)
        self.couplings = Couplings.from_dict(update_theory(self.theory))
        self.path_to_F2 = fiatlux_runcard["path_to_F2"]
        self.path_to_FL = fiatlux_runcard["path_to_FL"]
    
    def F2LO(self, x, Q):
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
        return self.couplings.a(q**2)[1] * 4 * np.pi

    def photon_fitting_scale(self, xgrid):
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
            return np.zeros(len(xgrid))
        
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

        photon_100GeV = np.array(
        [lux.EvaluatePhoton(x, self.q_in2).total / x for x in xgrid]
        )
        # TODO: fiatlux returns gamma(x) or x*gamma(x) ?
        
        pdfs = np.zeros((len(output.rotations.inputpids), len(xgrid)))
        for j, pid in enumerate(output.rotations.inputpids):
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
        for q2, elem in output.items():
            pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
            # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

        photon_fitting_scale = pdf_final[ph_id]

        # we want x * gamma(x)
        return xgrid * photon_fitting_scale
        