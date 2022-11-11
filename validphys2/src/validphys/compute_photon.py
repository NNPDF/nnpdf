"""Script that calls fiatlux to add the photon PDF."""

from validphys import lhapdfset
from validphys.api import API
import fiatlux
import numpy as np
from eko.couplings import Couplings
from eko.compatibility import update_theory
from eko.runner import Runner

import yaml

def photon_1GeV(xgrid, theoryid, fiatlux_runcard):
    r"""
    Compute the photon PDF for every point in the grid xgrid.

    Parameters
    ----------
        xgrid: list
            grid of the x points
        pdf_name: string
            name of the QCD set
        replica: int
            number of replica
    
    Returns
    -------
        photon_1GeV: numpay.array
            photon PDF at the scale 1 GeV
    """
    if fiatlux_runcard is None :
        return np.zeros(len(xgrid))
    set = lhapdfset.LHAPDFSet(fiatlux_runcard["pdf_name"], "replicas")
    qcd_pdfs = set.members[0] #use for the moment central replica TODO: change it

    def f2(x, q):
        # call yadism to give f2
        return 0.

    def fl(x, q):
        # call yadism to give fl
        return 0.

    def f2lo(x, q):
        # call yadism to give f20
        return 0.

    theory = API.theoryid(theoryid = theoryid).get_description().copy()
    # import ipdb; ipdb.set_trace()
    theory["nfref"] = None
    theory["nf0"] = None
    theory["fact_to_ren_scale_ratio"] = 1.
    theory["ModSV"] = None
    theory["IC"]=0
    theory["IB"]=0
    theory["FNS"] = "VFNS"
    q_in = 100
    q_in2 = q_in ** 2
    q_fin = theory["Q0"]
    theory["Q0"]= q_in

    lux = fiatlux.FiatLux(fiatlux_runcard)
    # TODO : we are passing a dict but fiatlux wants a yaml file
    qref_qed = theory["Qref"]

    couplings = Couplings.from_dict(update_theory(theory))
    def alpha_em(q):
        return couplings.a(q**2)[1] * 4 * np.pi
    
    lux.PlugAlphaQED(alpha_em, qref_qed)
    lux.PlugStructureFunctions(f2, fl, f2lo)
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

    photon_100GeV = np.array(
      [lux.EvaluatePhoton(x, q2).total / x for x in xgrid]
    )
    # TODO: fiatlux returns gamma(x) or x*gamma(x) ?
    runner = Runner(theory, operator_card)
    output = runner.get_output()
    pdfs = np.zeros((len(output.rotations.inputpids), len(xgrid)))
    for j, pid in enumerate(output.rotations.inputpids):
        if pid == 22 :
            pdfs[j] = photon_100GeV
        if not qcd_pdfs.hasFlavor(pid):
            continue
        pdfs[j] = np.array(
            [
                qcd_pdfs.xfxQ2(pid, x, q_in2) / x
                for x in xgrid
            ]
        )
    
    pdfs[11] = photon_100GeV
    for q2, elem in output.items():
        pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
        # error_final = np.einsum("ajbk,bk", elem.error, pdfs)

    photon_1GeV = pdf_final[11]

    return photon_1GeV