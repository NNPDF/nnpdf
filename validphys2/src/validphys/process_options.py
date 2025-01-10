"""
    Module to hold process dependent options

    Only variables included in the `_Vars` enum and processes included in the ``Processes`` dictionary are allowed.
"""

import dataclasses
from typing import Callable, Optional, Union

import numpy as np
from validobj.custom import Parser


class _Vars:
    x = "x"
    Q2 = "Q2"
    Q = "Q"
    y = "y"
    pT = "pT"
    ET = "ET"
    sqrts = "sqrts"
    ystar = "ystar"
    ydiff = "ydiff"
    m_jj = "m_jj"
    pT2 = "pT2"
    y_t = "y_t"
    y_ttBar = "y_ttBar"
    m_t2 = "m_t2"
    pT_t = "pT_t"
    m_ttBar = "m_ttBar"
    eta = "eta"
    abs_eta = "abs_eta"
    m_W2 = "m_W2"
    m_Z2 = "m_Z2"
    m_V2 = "m_V2"
    abs_eta_1 = "abs_eta_1"
    abs_eta_2 = "abs_eta_2"
    eta_1 = "eta_1"
    eta_2 = "eta_2"
    m_ll = "m_ll"
    m_ll2 = "m_ll2"
    abs_y = "abs_y"


class _KinematicsInformation:
    """Read the 3 columns dataframe corresponding to the values set in the
    ``kinematic_coverage`` field into a dictionary defining the name of the variables.

    Adds the special "sqrts" key unless it is already part of the kinematic coverage.

    Provides a ``.get_one_of`` method that accepts any number of variables
    """

    def __init__(self, kin_df, metadata):
        kin_cov = metadata.kinematic_coverage
        kins = {}
        for kin_lab, kin_values in zip(kin_cov, kin_df.values.T):
            kins[kin_lab] = kin_values

        if "sqrts" not in kin_cov:
            kins[_Vars.sqrts] = metadata.cm_energy

        self._kins = kins
        self._variables = list(kins.keys())

    def get_one_of(self, *variables):
        """Accepts any number of variables, returns one of them
        This is used for processes which might use slightly different definitions
        for instance ``pT`` vs ``ET`` but that can be used in the same manner"""
        for var in variables:
            if var in self._kins:
                return self._kins[var]
        raise KeyError(f"Need one of the following variables {self._variables} to continue")

    def __getitem__(self, key):
        return self.get_one_of(key)


@dataclasses.dataclass(frozen=True)
class _Process:
    name: str
    description: str
    accepted_variables: tuple[str]
    xq2map_function: Optional[Callable] = None

    def __hash__(self):
        return hash(self.name)

    def are_accepted_variables(self, kin_cov):
        """Check if the kinematic variables from the kinematic coverage are the same
        of the accepted variables."""
        # Accepting in any case the legacy variables
        if kin_cov == ["k1", "k2", "k3"]:
            return True
        # We check if kin_cov is a subset of self.accepted_variables
        kin_cov = [v for v in kin_cov if not v.startswith("extra_")]
        return set(self.accepted_variables).union(set(kin_cov)) == set(self.accepted_variables)

    def xq2map(self, kin_df, metadata):
        """Transform the kinematics dataframe into a x q dataframe
        For double hadronic processes the number of output point will be 2x ninput
        These functions should have access to both the kinematic dataframe and the
        metadata of the commondata
        """
        # Remove ``extra labels`` from kin_df
        if metadata.plotting.extra_labels is not None:
            for extra_label in metadata.plotting.extra_labels:
                kin_df = kin_df.drop(columns=extra_label)

        # Check if the kinematic variables defined in metadata corresponds to the
        # accepted variables
        if not self.are_accepted_variables(metadata.kinematic_coverage):
            raise NotImplementedError(
                f"kinematic variables are not supported for process {self.name}. You are using {metadata.kinematic_coverage}, please use {self.accepted_variables} ({metadata.name})"
            )

        if self.xq2map_function is None:
            raise NotImplementedError(f"xq2map is not implemented for {self.name}")

        # check that all variables in the dataframe are accepted by this process
        kininfo = _KinematicsInformation(kin_df, metadata)
        try:
            return self.xq2map_function(kininfo)
        except KeyError as e:
            raise NotImplementedError(
                f"Error trying to compute xq2map for process {self.name} ({metadata.name})"
            ) from e

    def __str__(self):
        return self.name


def _dis_xq2map(kin_info):
    """Variables in the dataframe should be x, Q2, y
    TODO: Once old variables are removed, remove if condition
    """
    x = kin_info.get_one_of("k1", _Vars.x)
    if "k2" in kin_info._kins:
        q2 = kin_info.get_one_of("k2") ** 2
    else:
        q2 = kin_info.get_one_of(_Vars.Q2)
    return x, q2


def _jets_xq2map(kin_info):
    # Then compute x, Q2
    pT = kin_info[_Vars.pT]
    ratio = pT / kin_info[_Vars.sqrts]
    rap = kin_info.get_one_of(_Vars.y, _Vars.eta, _Vars.abs_eta)
    x1 = 2 * ratio * np.exp(rap)
    x2 = 2 * ratio * np.exp(-rap)
    q2 = pT * pT
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))


def _shp_xq2map(kin_info):
    # Then compute x, Q2
    pT = kin_info[_Vars.pT]
    ratio = pT / kin_info[_Vars.sqrts]
    rap = kin_info[_Vars.eta]
    x1 = 2 * ratio * np.exp(rap)
    x2 = 2 * ratio * np.exp(-rap)
    q2 = pT * pT
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))


def _pht_xq2map(kin_info):
    # Then compute x, Q2
    ET = kin_info.get_one_of(_Vars.ET)
    ET2 = ET**2
    eta = kin_info.get_one_of(_Vars.eta, _Vars.y)
    sqrts = kin_info[_Vars.sqrts]

    # eta = y for massless particles
    x1 = ET / sqrts * np.exp(-eta)
    x2 = ET / sqrts * np.exp(eta)
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((ET2, ET2))


def _dijets_xq2map(kin_info):
    # Here we can have either ystar or ydiff, but in either case we need to do the same
    ylab_1 = kin_info.get_one_of(_Vars.ystar, _Vars.ydiff, _Vars.eta_1, _Vars.abs_eta_1)
    ylab_2 = kin_info.get_one_of(_Vars.ystar, _Vars.ydiff, _Vars.eta_2, _Vars.abs_eta_2)
    # Then compute x, Q2
    ratio = kin_info[_Vars.m_jj] / kin_info[_Vars.sqrts]
    x1 = ratio * np.exp(ylab_1)
    x2 = ratio * np.exp(-ylab_2)
    q2 = kin_info[_Vars.m_jj] * kin_info[_Vars.m_jj]
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))


def _hqp_yq_xq2map(kin_info):
    # Compute x, Q2
    #
    # Theory predictions computed with HT/4 ~ mt/2 for rapidity distr.
    # see section 3 from 1906.06535
    # HT defined in Eqn. (1) of 1611.08609
    rapidity = kin_info.get_one_of(_Vars.y_t, _Vars.y_ttBar, "k1")
    q2 = kin_info.get_one_of(_Vars.m_t2, "k2")
    ratio = np.sqrt(q2) / kin_info[_Vars.sqrts]
    x1 = ratio * np.exp(rapidity)
    x2 = ratio * np.exp(-rapidity)
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2)) / 4


def _hqp_ptq_xq2map(kin_info):
    # Compute x, Q2
    #
    # At LO pt ~ ptb
    # ht = 2.*sqrt(m_t2 + pT_t2)
    Q = (kin_info[_Vars.m_t2] + kin_info[_Vars.pT_t] * kin_info[_Vars.pT_t]) ** 0.5 / 2
    return Q / kin_info[_Vars.sqrts], Q * Q


def _hqp_mqq_xq2map(kin_info):
    # Compute x, Q2
    #
    # Theory predictions computed with HT/4 ~ m_ttbar/4
    Q = kin_info[_Vars.m_ttBar] / 4
    return Q / kin_info[_Vars.sqrts], Q * Q


def _inc_xq2map(kin_info):
    # Compute x, Q2
    # k2 necessary to take the mass for DY inclusive cross sections still not migrated
    mass2 = kin_info.get_one_of(_Vars.m_W2, _Vars.m_Z2, _Vars.m_t2, "k2")
    return np.sqrt(mass2) / kin_info[_Vars.sqrts], mass2


def _displusjet_xq2map(kin_info):
    """Computes x and q2 mapping for a DIS + J (J) process
    Uses Q2 as provided by the dictionary of kinematics variables
    and x = Q**4 / s / (pt**2 - Q**2)
    """
    q2 = kin_info[_Vars.Q2]
    # Consider ET and pT as equivalent for the purposes of the xq2 plot
    pt = kin_info.get_one_of(_Vars.ET, _Vars.pT)
    s = kin_info[_Vars.sqrts] ** 2
    x = q2 * q2 / s / (pt**2 - q2)
    return x, q2


def _dyboson_xq2map(kin_info):
    """
    Computes x and q2 mapping for pseudo rapidity observables
    originating from a W boson DY process.
    """
    mass2 = kin_info.get_one_of(_Vars.m_W2, _Vars.m_Z2, _Vars.m_V2)
    eta = kin_info.get_one_of(_Vars.eta, _Vars.y, _Vars.abs_eta)
    sqrts = kin_info[_Vars.sqrts]

    # eta = y for massless particles
    x1 = np.sqrt(mass2) / sqrts * np.exp(-eta)
    x2 = np.sqrt(mass2) / sqrts * np.exp(eta)
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((mass2, mass2))


def _dybosonpt_xq2map(kin_dict):
    """Compute x and q2 mapping for DY Z or W -> 2 leptons + jet process.
    Here pT refers to the transverse momentum of the boson.
    """
    pT = kin_dict[_Vars.pT]
    m_Z2 = kin_dict.get_one_of(_Vars.m_Z2, _Vars.m_W2, _Vars.m_ll2)

    sqrts = kin_dict[_Vars.sqrts]
    ET2 = m_Z2 + pT * pT
    x = (np.sqrt(ET2) + pT) / sqrts
    return x, ET2


def _dybosonptrap_xq2map(kin_info):
    """
    Computes x and q2 mapping for DY Z or W -> 2 leptons + jet process
    using the rapidity of the final lepton pair.
    """
    pT = kin_info[_Vars.pT]
    eta = kin_info.get_one_of(_Vars.eta, _Vars.y, _Vars.abs_y)
    m_ll2 = kin_info.get_one_of(_Vars.m_ll2, _Vars.m_Z2)
    sqrts = kin_info[_Vars.sqrts]
    ET2 = m_ll2 + pT * pT
    x1 = (np.sqrt(ET2) + pT) / sqrts * np.exp(-eta)
    x2 = (np.sqrt(ET2) + pT) / sqrts * np.exp(+eta)
    x = np.concatenate((x1, x2))
    return x, np.concatenate((ET2, ET2))


def _singletop_xq2map(kin_dict):

    y_t = kin_dict[_Vars.y_t]
    sqrts = kin_dict[_Vars.sqrts]
    m_t2 = kin_dict[_Vars.m_t2]

    q2 = m_t2
    ratio = np.sqrt(q2) / sqrts
    x1 = ratio * np.exp(y_t)
    x2 = ratio * np.exp(-y_t)
    x = np.concatenate((x1, x2))
    return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))


def _dymll_xq2map(kin_info):
    """
    Computes x and q2 mapping for DY Z -> 2 leptons mass.
    x is approximated as x = sqrt(x1*x2) with m_ll^2 = x1*x2*s
    """

    m_ll = kin_info.get_one_of(_Vars.m_ll)
    sqrts = kin_info.get_one_of(_Vars.sqrts)
    m_ll2 = m_ll**2
    x = m_ll / sqrts

    return x, m_ll2


DIS = _Process(
    "DIS",
    "Deep Inelastic Scattering",
    accepted_variables=(_Vars.x, _Vars.Q2, _Vars.y, _Vars.Q),
    xq2map_function=_dis_xq2map,
)

JET = _Process(
    "JET",
    "Single Jet production",
    accepted_variables=(_Vars.y, _Vars.pT, _Vars.sqrts, _Vars.pT2),
    xq2map_function=_jets_xq2map,
)

SHP = _Process(
    "SHP",
    "Single Hadron Production",
    accepted_variables=(_Vars.eta, _Vars.pT, _Vars.sqrts),
    xq2map_function=_shp_xq2map,
)

DIJET = _Process(
    "DIJET",
    "DiJets production",
    accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff),
    xq2map_function=_dijets_xq2map,
)

JET_POL = _Process(
    "JET_POL",
    "Longitudinal double-spin asymmetry in inclusive jet production",
    accepted_variables=(_Vars.eta, _Vars.pT, _Vars.sqrts, _Vars.abs_eta),
    xq2map_function=_jets_xq2map,
)

DIJET_POL = _Process(
    "DIJET_POL",
    "Longitudinal double-spin asymmetry in dijets production",
    accepted_variables=(
        _Vars.m_jj,
        _Vars.sqrts,
        _Vars.abs_eta_2,
        _Vars.abs_eta_1,
        _Vars.eta_1,
        _Vars.eta_2,
    ),
    xq2map_function=_dijets_xq2map,
)

HQP_YQ = _Process(
    "HQP_YQ",
    "(absolute) rapidity of top quark in top pair production",
    accepted_variables=(
        _Vars.y_t,
        _Vars.y_ttBar,
        _Vars.m_t2,
        _Vars.sqrts,
        _Vars.m_ttBar,
        _Vars.pT_t,
    ),
    xq2map_function=_hqp_yq_xq2map,
)

HQP_PTQ = _Process(
    "HQP_PTQ",
    "Transverse momentum of top quark in top pair production",
    accepted_variables=(
        _Vars.pT_t,
        _Vars.m_ttBar,
        _Vars.y_t,
        _Vars.y_ttBar,
        _Vars.sqrts,
        _Vars.m_t2,
    ),
    xq2map_function=_hqp_ptq_xq2map,
)

HQP_MQQ = _Process(
    "HQP_MQQ",
    "Invariant mass of top quark pair in top pair production",
    accepted_variables=(_Vars.m_ttBar, _Vars.y_t, _Vars.y_ttBar, _Vars.sqrts, _Vars.m_t2),
    xq2map_function=_hqp_mqq_xq2map,
)

INC = _Process(
    "INC",
    "Inclusive cross section",
    accepted_variables=("zero", _Vars.sqrts, _Vars.m_W2, _Vars.m_Z2, _Vars.m_t2),
    xq2map_function=_inc_xq2map,
)

HERAJET = _Process(
    "HERAJET",
    "DIS + j production",
    accepted_variables=(_Vars.pT, _Vars.Q2, _Vars.sqrts, _Vars.ET),
    xq2map_function=_displusjet_xq2map,
)


DY_2L = _Process(
    "DY_2L",
    "DY W or Z -> 2 leptons ",
    accepted_variables=(
        _Vars.y,
        _Vars.eta,
        _Vars.m_W2,
        _Vars.m_Z2,
        _Vars.m_V2,
        _Vars.sqrts,
        _Vars.abs_eta,
    ),
    xq2map_function=_dyboson_xq2map,
)

DY_MLL = _Process(
    "DY_MLL",
    "DY Z -> ll mass of lepton pair",
    accepted_variables=(_Vars.m_ll, _Vars.sqrts),
    xq2map_function=_dymll_xq2map,
)

DY_PT = _Process(
    "DY_PT",
    "DY W or Z (2 leptons) + j boson transverse momentum",
    accepted_variables=(_Vars.pT, _Vars.m_W2, _Vars.m_Z2, _Vars.sqrts, _Vars.y, _Vars.m_ll2),
    xq2map_function=_dybosonpt_xq2map,
)

DY_PT_RAP = _Process(
    "DY_PT",
    "DY W or Z (2 leptons) + j boson transverse momentum",
    accepted_variables=(
        _Vars.pT,
        _Vars.m_W2,
        _Vars.m_Z2,
        _Vars.sqrts,
        _Vars.y,
        _Vars.abs_y,
        _Vars.eta,
        _Vars.m_ll2,
    ),
    xq2map_function=_dybosonptrap_xq2map,
)


POS_XPDF = _Process("POS_XPDF", "Positivity of MS bar PDFs", accepted_variables=(_Vars.x, _Vars.Q2))

POS_DIS = _Process(
    "POS_DIS", "Positivity of F2 structure functions", accepted_variables=(_Vars.x, _Vars.Q2)
)

PHT = _Process(
    "PHT",
    "Photon production",
    accepted_variables=(_Vars.eta, _Vars.ET, _Vars.sqrts),
    xq2map_function=_pht_xq2map,
)

SINGLETOP = _Process(
    "SINGLETOP",
    "Single top production",
    accepted_variables=(_Vars.m_t2, _Vars.sqrts, _Vars.y_t, _Vars.pT_t),
    xq2map_function=_singletop_xq2map,
)


PROCESSES = {
    "DIS": DIS,
    "DIS_NC": dataclasses.replace(DIS, name="DIS_NC"),
    "DIS_CC": dataclasses.replace(DIS, name="DIS_CC"),
    "DIS_NCE": dataclasses.replace(DIS, name="DIS_NCE"),
    "DIS_POL": dataclasses.replace(DIS, name="DIS_POL"),
    "DIS_NC_CHARM": dataclasses.replace(DIS, name="DIS_NC_CHARM"),
    "DIS_NC_BOTTOM": dataclasses.replace(DIS, name="DIS_NC_BOTTOM"),
    "JET": JET,
    "DIJET": DIJET,
    "SHP_ASY": SHP,
    "HQP_YQ": HQP_YQ,
    "HQP_YQQ": dataclasses.replace(HQP_YQ, name="HQP_YQQ"),
    "HQP_PTQ": HQP_PTQ,
    "HQP_MQQ": HQP_MQQ,
    "INC": INC,
    "HERAJET": HERAJET,
    "HERADIJET": dataclasses.replace(HERAJET, name="HERADIJET", description="DIS + jj production"),
    "JET_POL": JET_POL,
    "DIJET_POL": DIJET_POL,
    "DY_Z_Y": dataclasses.replace(DY_2L, name="DY_Z_Y", description="DY Z -> ll (pseudo)rapidity"),
    "DY_MLL": DY_MLL,
    "DY_W_ETA": dataclasses.replace(
        DY_2L, name="DY_W_ETA", description="DY W -> l nu pseudorapidity"
    ),
    "DY_VB_ETA": dataclasses.replace(
        DY_2L, name="DY_VB_ETA", description="DY Z/W -> ll pseudorapidity"
    ),
    "DY_NC_PT": dataclasses.replace(DY_PT, name="DY_NC_PT", description="DY Z (ll) + j"),
    "DY_CC_PT": dataclasses.replace(DY_PT, name="DY_CC_PT", description="DY W + j"),
    "DY_NC_PTRAP": dataclasses.replace(DY_PT_RAP, name="DY_NC_PTRAP", description="DY Z (ll) + j"),
    "POS_XPDF": POS_XPDF,
    "POS_DIS": POS_DIS,
    "PHT": PHT,
    "SINGLETOP": SINGLETOP,
}


@Parser
def ValidProcess(process_name) -> Union[_Process, str]:
    return PROCESSES.get(process_name.upper(), process_name.upper())
