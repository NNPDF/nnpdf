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

def photon_fitting_scale(xgrid, theoryid, fiatlux_runcard, replica):
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
    if fiatlux_runcard is None :
        return np.zeros(len(xgrid))
    pdf_name = fiatlux_runcard["pdf_name"]
    qcd_pdfs = lhapdf.mkPDF(pdf_name, replica)

    theory = API.theoryid(theoryid = theoryid).get_description().copy()
    # theory["nfref"] = None
    # theory["nf0"] = None
    # theory["fact_to_ren_scale_ratio"] = 1.
    theory["ModSV"] = None
    # theory["IC"]=0
    # theory["IB"]=0
    theory["FNS"] = "VFNS"
    q_in = 100
    q_in2 = q_in ** 2
    q_fin = theory["Q0"]
    theory["Q0"]= q_in
    qref = theory["Qref"]

    # lux = fiatlux.FiatLux(fiatlux_runcard)
    # we have a dict but fiatlux wants a yaml file
    # TODO : remove this trick
    ff = open('fiatlux_runcard.yml', 'w+')
    yaml.dump(fiatlux_runcard, ff)

    lux = fiatlux.FiatLux('fiatlux_runcard.yml')

    couplings = Couplings.from_dict(update_theory(theory))
    def alpha_em(q):
        return couplings.a(q**2)[1] * 4 * np.pi
    
    lux.PlugAlphaQED(alpha_em, qref)
    
    path_to_F2 = fiatlux_runcard["path_to_F2"]
    path_to_FL = fiatlux_runcard["path_to_FL"]
    # xir = theory["XIR"]
    # xif = theory["XIF"]
    
    def F2LO(x, Q):
        mcharm = theory["mc"]
        mbottom = theory["mb"]
        mtop = theory["mt"]
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
        e2q = [e2d, e2u, e2d, e2u, e2d, e2u]
        res = 0
        for i in range(1, hq+1):
            res += e2q[i-1] * (qcd_pdfs.xfxQ(x, Q)[i] + qcd_pdfs.xfxQ(x, Q)[-i])
        return res

    f2 = StructureFunction(path_to_F2, qcd_pdfs)
    fl = StructureFunction(path_to_FL, qcd_pdfs)
    
    lux.PlugStructureFunctions(f2.FxQ, fl.FxQ, F2LO)
    
    lux.InsertInelasticSplitQ([4.18, 1e100])

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
                Q2grid=[q_fin**2],
                inputgrid=None,
                targetgrid=None,
                inputpids=None,
                targetpids=None,
            ).items()
        )
    )
    for x in xgrid:
        print(lux.EvaluatePhoton(x, q_in2).total / x)

    photon_100GeV = np.array(
      [lux.EvaluatePhoton(x, q_in2).total / x for x in xgrid]
    )
    # TODO: fiatlux returns gamma(x) or x*gamma(x) ?
    
    pdfs = np.zeros((len(output.rotations.inputpids), len(xgrid)))
    for j, pid in enumerate(output.rotations.inputpids):
        if pid == 22 :
            pdfs[j] = photon_100GeV
            ph_id = j
        if not qcd_pdfs.hasFlavor(pid):
            continue
        pdfs[j] = np.array(
            [
                qcd_pdfs.xfxQ2(pid, x, q_in2) / x
                for x in xgrid
            ]
        )

    runner = Runner(theory, operator_card)
    output = runner.get_output()
    for q2, elem in output.items():
        pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
        # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

    photon_fitting_scale = pdf_final[ph_id]

    # we want x * gamma(x)
    return xgrid * photon_fitting_scale