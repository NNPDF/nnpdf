"""
    Module to hold process dependent options

    Only variables included in the `_Vars` enum and processes included in the ``Processes`` enum are allowed
"""
import dataclasses
from typing import Callable, Optional, Tuple
import numpy as np

from validobj.custom import Parser

TMASS = 173.3

class _Vars:
    x = "x"
    Q2 = "Q2"
    Q = "Q"
    y = "y"
    pT = "pT" 
    sqrts = "sqrts"
    ystar = "ystar"
    ydiff = "ydiff"
    m_jj = "m_jj"
    p_T2 = "p_T2"
    y_t = "y_t"
    y_ttBar = "y_ttBar"
    m_t2 = "m_t2"
    pT_t = "pT_t"
    m_ttBar = "m_ttBar"

def map_to_metadata(kin_df, kin_cov):
    kins = {}
    for kin_lab, kin_values in zip(kin_cov, kin_df.values.T):
        kins[kin_lab] = kin_values
    return kins

@dataclasses.dataclass(frozen=True)
class _Process:
    name: str
    description: str
    accepted_variables: Tuple[str]
    xq2map_function: Optional[Callable] = None

    def __hash__(self):
        return hash(self.name)

    def same_kin_variables(self, kin_cov):
        """Check if the kinematic variables from the kinematic coverage are the same
        of the accepted variables."""
        # Accepting in any case the legacy variables
        if kin_cov == ["k1", "k2", "k3"]:
            return True
        # We check if kin_cov is a subset of self.accepted_variables
        return set(self.accepted_variables).union(set(kin_cov)) == set(self.accepted_variables)
    
    def xq2map(self, kin_df, metadata):
        """Transform the kinematics dataframe into a x q dataframe
        For double hadronic processes the number of output point will be 2x ninput
        These functions should have access to both the kinematic dataframe and the
        metadata of the commondata
        """
        # Remove extra labels from kin_df
        if metadata.plotting.extra_labels is not None:
            for extra_label in metadata.plotting.extra_labels:
                kin_df = kin_df.drop(columns=extra_label) 
        # Check if the kinematic variables defined in metadata corresponds to the 
        # accepted variables
        if not self.same_kin_variables(metadata.kinematic_coverage):
            raise NotImplementedError(f"kinematic variables are not supported for this process. You are using {metadata.kinematic_coverage}, please use {self.accepted_variables}")
        if self.xq2map_function is None:
            raise NotImplementedError(f"xq2map is not implemented for {self.name}")
        # check that all variables in the dataframe are accepted by this process
        return self.xq2map_function(map_to_metadata(kin_df, metadata.kinematic_coverage))

    def __str__(self):
        return self.name


def _dis_xq2map(kin_dict):
    """In the old style commondata, the variables in the dataframe were ``x, Q2, y``
    but due to the transformations that happen inside validphys they become ``x, Q, y``
    """
    x = kin_dict["k1"]
    q = kin_dict["k2"]
    return x, q * q
        
def _jets_xq2map(kin_dict):
    # Then compute x, Q2
    if "p_T2" in kin_dict:
        pT = np.sqrt(kin_dict["p_T2"])
    else:
        pT = kin_dict["pT"]
    ratio = pT / kin_dict["sqrts"]
    x1 = 2 * ratio * np.exp(kin_dict["y"])
    x2 = 2 * ratio * np.exp(-kin_dict["y"])
    q2 = pT * pT
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))

def _dijets_xq2map(kin_dict):
    # Here we can have either ystar or ydiff, but in either case we need to do the same
    ylab = "ystar" if "ystar" in kin_dict else "ydiff"
    # Then compute x, Q2
    ratio = kin_dict["m_jj"] / kin_dict["sqrts"]
    x1 = ratio * np.exp(kin_dict[ylab])
    x2 = ratio * np.exp(-kin_dict[ylab])
    q2 = kin_dict["m_jj"] * kin_dict["m_jj"]
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))

def _hqp_yq_xq2map(kin_dict):
    # Compute x, Q2
    if list(kin_dict.keys()) == ["k1", "k2", "k3"]:
        kin_dict["y_t"] = kin_dict["k1"]
        #k2 is just the mass of the top, not squared
        kin_dict["m_t2"] = kin_dict["k2"]**2
        kin_dict["sqrts"] = kin_dict["k3"]
    m_label = "m_t2" if "m_t2" in kin_dict else "m_ttBar"
    mass2 = kin_dict[m_label]**2 if m_label == "m_ttBar" else kin_dict[m_label]
    ratio = np.sqrt(mass2) / kin_dict["sqrts"]
    x1 = ratio * np.exp(kin_dict["y_t"])
    x2 = ratio * np.exp(-kin_dict["y_t"])
    q2 = mass2
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))

def _hqp_yqq_xq2map(kin_dict):
    # Compute x, Q2
    ratio = np.sqrt(kin_dict["m_t2"]) / kin_dict["sqrts"]
    x1 = ratio * np.exp(kin_dict["y_ttBar"])
    x2 = ratio * np.exp(-kin_dict["y_ttBar"])
    q2 = kin_dict["m_t2"]
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))

def _hqp_ptq_xq2map(kin_dict):
    # Compute x, Q2
    QMASS2 = TMASS * TMASS
    Q = np.sqrt(QMASS2 + kin_dict["pT_t"] * kin_dict["pT_t"]) + kin_dict["pT_t"]
    return Q / kin_dict["sqrts"], Q * Q

DIS = _Process(
    "DIS",
    "Deep Inelastic Scattering",
    accepted_variables=(_Vars.x, _Vars.Q2, _Vars.y, _Vars.Q),
    xq2map_function=_dis_xq2map,
)

JET = _Process(
    "JET",
    "Single Jet production",
    accepted_variables=(_Vars.y, _Vars.pT, _Vars.sqrts, _Vars.p_T2),
    xq2map_function=_jets_xq2map,
    
)

DIJET = _Process(
    "DIJET",
    "DiJets Production",
    accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff),
    xq2map_function=_dijets_xq2map,   
)

HQP_YQ = _Process(
    "HQP_YQ",
    "Normalized differential cross section w.r.t. absolute rapidity of t",
    accepted_variables=(_Vars.y_t, _Vars.m_t2, _Vars.sqrts, _Vars.m_ttBar),
    xq2map_function=_hqp_yq_xq2map
)

HQP_YQQ = _Process(
    "HQP_YQQ",
    "Differential cross section w.r.t. absolute rapidity of ttBar",
    accepted_variables=(_Vars.y_ttBar, _Vars.m_t2, _Vars.sqrts),
    xq2map_function=_hqp_yqq_xq2map
)

HQP_PTQ = _Process(
    "HQP_PTQ",
    "Normalized double differential cross section w.r.t. absolute rapidity and transverse momentum of t",
    accepted_variables=(_Vars.pT_t, _Vars.y_t, _Vars.sqrts),
    xq2map_function=_hqp_ptq_xq2map
)

PROCESSES = {
    "DIS": DIS,
    "DIS_NC": dataclasses.replace(DIS, name="DIS_NC"),
    "DIS_CC": dataclasses.replace(DIS, name="DIS_CC"),
    "DIS_NCE": dataclasses.replace(DIS, name="DIS_NCE"),
    "JET": JET,
    "DIJET": DIJET,
    "HQP_YQ": HQP_YQ,
    "HQP_YQQ": HQP_YQQ,
    "HQP_PTQ": HQP_PTQ,
}


@Parser
def ValidProcess(process_name) -> _Process:
    return PROCESSES.get(process_name.upper(), process_name.upper())
