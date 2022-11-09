"""Script that calls fiatlux to add the photon PDF."""

from validphys import lhapdfset
from validphys.api import API
import fiatlux
import numpy as np
from eko.couplings import Couplings
from eko.compatibility import update_theory
from eko.runner import Runner

import yaml

def photon_1GeV(xgrid, theoryID, pdf_name, replica):
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
    set = lhapdfset.LHAPDFSet(pdf_name, "replicas")
    pdfs = np.zeros((len(set.flavors()), len(xgrid)))
    for j, pid in enumerate(set.flavors()):
        if not set.hasFlavor(pid):
            continue
        pdfs[j] = np.array(
            [
                set.xfxQ(x, 100., replica, pid) / x
                for x in xgrid
            ]
        )
    out_grid = {}

    def f2(x, q):
        # call yadism to give f2
        return 0.

    def fl(x, q):
        # call yadism to give fl
        return 0.

    def f2lo(x, q):
        # call yadism to give f20
        return 0.

    theory = API.theoryid(theoryid = theoryID).get_description()
    # fiatlux_runcard = <load runcard>
    fiatlux_runcard = {
        'apfel': True,
        'qed_running': True,
        'q2_max': '1e9',
        'eps_base': '1e-5',
        'eps_rel': '1e-1',
        'mproton': 0.938272046,
        'mum_proton': 2.792847356,
        'elastic_param': 'A1_world_pol_spline',
        'elastic_electric_rescale': 1,
        'elastic_magnetic_rescale': 1,
        'inelastic_param': 'LHAPDF_Hermes_ALLM_CLAS',
        'rescale_r_twist4': 0,
        'rescale_r': 1, 'allm_limits': 0,
        'rescale_non_resonance': 1,
        'rescale_resonance': 1,
        'use_mu2_as_upper_limit': False,
        'q2min_inel_override': 0.0,
        'q2max_inel_override': '1E300',
        'lhapdf_transition_q2': 9,
        'verbose': True
    } 
    # TODO : move this dict into the runcard and load it
    lux = fiatlux.FiatLux(yaml.dump(fiatlux_runcard))
    qref_qed = theory["Qref"]

    couplings = Couplings.from_dict(update_theory(theory))
    def alpha_em(q):
        return couplings.a(q**2)[1] * 4 * np.pi
    
    lux.PlugAlphaQED(alpha_em, qref_qed)
    lux.PlugStructureFunctions(f2, fl, f2lo)
    lux.InsertInelasticSplitQ([4.18, 1e100])
    eko_operator = dict(
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
                Q2grid=[1.**2],
                inputgrid=None,
                targetgrid=None,
                inputpids=None,
                targetpids=None,
            ).items()
        )
    )
    q2 = 100.**2
    photon_100GeV = np.array([])
    for x in xgrid:
        pht = lux.EvaluatePhoton(x, q2)
        np.append(photon_100GeV, pht.total)
    pdfs[11] = photon_100GeV

    runner = Runner(theory, eko_operator)
    eko = runner.get_output()

    for q2, elem in eko.items():
        pdf_final = np.einsum("ajbk,bk", elem.operator, pdfs)
        error_final = np.einsum("ajbk,bk", elem.error, pdfs)
        out_grid[q2] = {
            "pdfs": dict(zip(eko.rotations.targetpids, pdf_final)),
            "errors": dict(zip(eko.rotations.targetpids, error_final)),
        }
    photon_1GeV = pdf_final[11]

    return photon_1GeV