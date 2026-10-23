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
    abs_y = "abs_y"
    pT = "pT"
    ET = "ET"
    sqrts = "sqrts"
    ystar = "ystar"
    ydiff = "ydiff"
    ymax = "ymax"
    yb = "yb"
    m_jj = "m_jj"
    pT2 = "pT2"
    pTavg = "pTavg"
    y_t = "y_t"
    y_ttBar = "y_ttBar"
    m_t2 = "m_t2"
    m_t = "m_t"
    pT_t = "pT_t"
    m_ttBar = "m_ttBar"
    eta = "eta"
    abs_eta = "abs_eta"
    m_W2 = "m_W2"
    m_Z2 = "m_Z2"
    m_V2 = "m_V2"
    m_W = "m_W"
    m_Z = "m_Z"
    M2 = "M2"
    abs_eta_1 = "abs_eta_1"
    abs_eta_2 = "abs_eta_2"
    eta_1 = "eta_1"
    eta_2 = "eta_2"
    m_ll = "m_ll"
    m_ll2 = "m_ll2"
    abs_y = "abs_y"
    # Variables used by fragmentation processes (SIA, SIDIS)
    z = "z"
    xp = "xp"
    ph = "ph"
    xi = "xi"


class _KinematicsInformation:
    """Read the 3 columns dataframe corresponding to the values set in the
    ``kinematic_coverage`` field into a dictionary defining the name of the variables.

    Adds the special "sqrts" key unless it is already part of the kinematic coverage.

    Provides a ``.get_one_of`` method that accepts any number of variables

    The process field of the dataset name is available as ``.dataset_process``, for the maps that
    need more information than the kinematic variables (e.g., the identity of the observed hadron).
    It is not a kinematic variable.
    """

    def __init__(self, kin_df, metadata):
        self.dataset_process = metadata.process
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
    zq2map_function: Optional[Callable] = None
    # Kinematic maps this process can be mapped to: "xq2", "zq2" or both.
    # If not given, it is inferred from the map functions that have been defined.
    kinematic_maps: Optional[tuple[str, ...]] = None

    def __post_init__(self):
        map_functions = {"xq2": self.xq2map_function, "zq2": self.zq2map_function}
        maps = self.kinematic_maps
        if maps is None:
            maps = tuple(name for name, function in map_functions.items() if function is not None)
        maps = tuple(maps)

        unknown = set(maps) - set(map_functions)
        if unknown:
            raise ValueError(
                f"Unknown kinematic maps {sorted(unknown)} for process {self.name}, "
                f"the options are {list(map_functions)}"
            )
        for name, function in map_functions.items():
            if (name in maps) != (function is not None):
                raise ValueError(
                    f"The kinematic maps of process {self.name} are {maps} "
                    f"which is inconsistent with the {name}map_function defined"
                )
        # the dataclass is frozen
        object.__setattr__(self, "kinematic_maps", maps)

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
        # sqrts is always accepted
        accepted = set(self.accepted_variables).union({_Vars.sqrts})
        return accepted.union(set(kin_cov)) == accepted

    def _apply_map(self, kin_df, metadata, map_function, map_name):
        """Validate the kinematics and apply ``map_function`` to them.

        Shared by all the kinematic maps (e.g., ``xq2map``). ``map_name`` is only used
        to build the error messages.
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

        if map_function is None:
            raise NotImplementedError(f"{map_name} is not implemented for {self.name}")

        # check that all variables in the dataframe are accepted by this process
        kininfo = _KinematicsInformation(kin_df, metadata)
        try:
            return map_function(kininfo)
        except KeyError as e:
            raise NotImplementedError(
                f"Error trying to compute {map_name} for process {self.name} ({metadata.name})"
            ) from e

    def xq2map(self, kin_df, metadata):
        """Transform the kinematics dataframe into a x q dataframe
        For double hadronic processes the number of output point will be 2x ninput
        These functions should have access to both the kinematic dataframe and the
        metadata of the commondata
        """
        return self._apply_map(kin_df, metadata, self.xq2map_function, "xq2map")

    def zq2map(self, kin_df, metadata):
        """Transform the kinematics dataframe into a z q dataframe
        This is the equivalent of ``xq2map`` for processes involving a fragmentation function,
        where z is the momentum fraction of the observed hadron.
        """
        return self._apply_map(kin_df, metadata, self.zq2map_function, "zq2map")

    def __str__(self):
        return self.name


def _dis_xq2map(kin_info):
    """Variables in the dataframe should be x, Q2, y
    TODO: Once old variables are removed, remove if condition
    """
    x = kin_info.get_one_of("k1", _Vars.x)
    if "k2" in kin_info._kins:
        q2 = kin_info.get_one_of("k2") ** 2
    elif _Vars.Q in kin_info._kins:
        q = kin_info.get_one_of(_Vars.Q)
        q2 = q**2
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
    # Here we can have either ystar or ymax or ydiff, but in either case we need to do the same
    if _Vars.m_jj in kin_info._kins:
        ylab_1 = kin_info.get_one_of(_Vars.ystar, _Vars.ydiff, _Vars.ymax, _Vars.eta_1, _Vars.abs_eta_1)
        ylab_2 = kin_info.get_one_of(_Vars.ystar, _Vars.ydiff, _Vars.ymax, _Vars.eta_2, _Vars.abs_eta_2)
        ratio = kin_info[_Vars.m_jj] / kin_info[_Vars.sqrts]
        x1 = ratio * np.exp(ylab_1)
        x2 = ratio * np.exp(-ylab_2)
        q2 = kin_info[_Vars.m_jj] * kin_info[_Vars.m_jj]
    elif _Vars.pTavg in kin_info._kins:
        m_jjlab = 2 * kin_info[_Vars.pTavg] * np.cosh(kin_info[_Vars.ystar])
        ratio = m_jjlab / kin_info[_Vars.sqrts]
        yblab = kin_info[_Vars.yb]
        x1 = ratio * np.exp(yblab)
        x2 = ratio * np.exp(-yblab)
        q2 = m_jjlab * m_jjlab
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
    try:
        mass = kin_info.get_one_of(_Vars.m_W, _Vars.m_Z, _Vars.m_t, _Vars.m_ll)
        mass2 = mass**2
    except KeyError:
        mass2 = kin_info.get_one_of(_Vars.m_W2, _Vars.m_Z2, _Vars.m_t2, _Vars.m_ll)
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
    try:
        mass2 = kin_info.get_one_of(
            _Vars.m_W2, _Vars.m_Z2, _Vars.m_V2, _Vars.m_ll2, _Vars.M2
        )
    except KeyError:
        mass = kin_info.get_one_of(_Vars.m_W, _Vars.m_Z, _Vars.m_ll, _Vars.m_Z)
        mass2 = mass**2
    eta = kin_info.get_one_of(_Vars.eta, _Vars.y, _Vars.abs_eta, _Vars.abs_y)
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

    try:
        m_Z2 = kin_dict.get_one_of(_Vars.m_Z2, _Vars.m_W2, _Vars.m_ll2)
    except KeyError:
        m_Z2 = kin_dict.get_one_of(_Vars.m_ll) ** 2

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
    try:
        m_ll2 = kin_info.get_one_of(_Vars.m_Z2, _Vars.m_W2, _Vars.m_ll2)
    except KeyError:
        m_ll2 = kin_info.get_one_of(_Vars.m_ll) ** 2
        
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


def _sidis_zq2map(kin_info):
    """Computes z and q2 mapping for a SIDIS process.
    z is the momentum fraction of the observed hadron and the scale is
    taken from either Q or Q2 (as in DIS).
    """
    z = kin_info[_Vars.z]
    if _Vars.Q in kin_info._kins:
        q2 = kin_info[_Vars.Q] ** 2
    else:
        q2 = kin_info[_Vars.Q2]
    return z, q2


# SIA datasets are deliberately separated by hadron type: the mass of the observed hadron enters
# the z-Q2 map when z is derived from a momentum. The hadron is read from the process field of the
# dataset name, which is expected to be of the form ``SIA-<HADRON>``, with <HADRON> one of the
# keys below. Masses in GeV, taken from the PDG.
# TODO (open question): should the masses instead be read from the theory card used for the
# fragmentation functions, so that the value used here is always the one used in the theory?
_HADRON_MASSES = {"PI": 0.139570, "KA": 0.493677, "PRO": 0.938272}


def _hadron_mass(kin_info):
    """Mass of the observed hadron, from the process field of the dataset (``SIA-<HADRON>``)"""
    hadron = kin_info.dataset_process.split("-")[-1]
    if hadron not in _HADRON_MASSES:
        raise KeyError(
            f"Unknown hadron {hadron!r} in process {kin_info.dataset_process!r}, "
            f"the options are {list(_HADRON_MASSES)}"
        )
    return _HADRON_MASSES[hadron]


def _sia_sqrts(kin_info):
    """Center of mass energy of an SIA observable, which sets the scale Q = sqrts"""
    sqrts = kin_info[_Vars.sqrts]
    if sqrts is None:
        raise KeyError("The center of mass energy is needed but it is not fixed for this dataset")
    return sqrts


def _sia_z_from_momentum(momentum, mass, sqrts):
    """z = 2 E_h / sqrts, with E_h the energy of a hadron of the given mass and momentum"""
    return 2 * np.sqrt(momentum**2 + mass**2) / sqrts


def _sia_zq2(z, sqrts):
    """Return z together with the scale Q2 = sqrts^2 (one value per point, like z)"""
    return z, sqrts**2 * np.ones_like(z)


# The different ways in which the momentum (or energy) fraction of the observed hadron can be given
_SIA_VARIABLES = (_Vars.z, _Vars.xp, _Vars.ph, _Vars.xi)


def _sia_zq2map(kin_info):
    """Computes z and q2 mapping for an SIA process.
    The observable is measured as a function of one (and only one) of the following variables:

    - z, the energy of the hadron scaled to the beam energy, z = 2 E_h / sqrts
    - xp, the momentum of the hadron scaled to the beam momentum, p = xp * sqrts / 2
    - ph, the momentum of the hadron (in GeV)
    - xi = ln(1/xp)

    The mass of the hadron is needed to obtain z from a momentum, i.e., in all but the first case.
    The scale is the center of mass energy, Q2 = sqrts^2.
    """
    given = [variable for variable in _SIA_VARIABLES if variable in kin_info._kins]
    if len(given) != 1:
        raise KeyError(f"Exactly one of {list(_SIA_VARIABLES)} is needed, the dataset has {given}")
    (variable,) = given

    sqrts = _sia_sqrts(kin_info)
    if variable == _Vars.z:
        return _sia_zq2(kin_info[_Vars.z], sqrts)
    if variable == _Vars.xp:
        momentum = kin_info[_Vars.xp] * sqrts / 2
    elif variable == _Vars.ph:
        momentum = kin_info[_Vars.ph]
    else:
        momentum = np.exp(-kin_info[_Vars.xi]) * sqrts / 2
    z = _sia_z_from_momentum(momentum, _hadron_mass(kin_info), sqrts)
    return _sia_zq2(z, sqrts)


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
    accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff, _Vars.ymax),
    xq2map_function=_dijets_xq2map,
)

DIJET_3D = _Process(
    "DIJET_3D",
    "DiJets production where the measured quantity is triple differential cross section",
    accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff, _Vars.ymax, _Vars.yb, _Vars.pTavg),
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
    accepted_variables=("zero", _Vars.sqrts, _Vars.m_W2, _Vars.m_Z2, _Vars.m_t2, _Vars.m_ll),
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
        _Vars.abs_y,
        _Vars.eta,
        _Vars.m_W2,
        _Vars.m_Z2,
        _Vars.m_V2,
        _Vars.m_ll,
        _Vars.m_ll2,
        _Vars.sqrts,
        _Vars.abs_eta,
        _Vars.m_ll2,
        _Vars.M2,
        _Vars.m_Z
    ),
    xq2map_function=_dyboson_xq2map,
)

DY_MLL = _Process(
    "DY_MLL",
    "DY Z -> ll mass of lepton pair",
    accepted_variables=(_Vars.m_ll, _Vars.sqrts, _Vars.y),
    xq2map_function=_dymll_xq2map,
)

DY_PT = _Process(
    "DY_PT",
    "DY W or Z (2 leptons) + j boson transverse momentum",
    accepted_variables=(_Vars.pT, _Vars.m_W2, _Vars.m_Z2, _Vars.sqrts, _Vars.y, _Vars.m_ll2, _Vars.m_ll),
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
        _Vars.m_ll,
    ),
    xq2map_function=_dybosonptrap_xq2map,
)


POS_XPDF = _Process("POS_XPDF", "Positivity of MS bar PDFs", accepted_variables=(_Vars.x, _Vars.Q2))

POS_DIS = _Process(
    "POS_DIS", "Positivity of F2 structure functions", accepted_variables=(_Vars.x, _Vars.Q2)
)

POS_POL = _Process(
    "POS_POL",
    "Polarized DIS positivity boundary condition",
    accepted_variables=(_Vars.x, _Vars.Q2, _Vars.y),
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

# TODO (open question): SIDIS is a single process here, with no split by current or by hadron.
# DIS is split into NC/CC variants in ``PROCESSES`` and the SIA datasets are deliberately
# separated by hadron type, since the hadron mass enters their z-Q2 map. Nothing needs that
# for SIDIS yet (the maps only use x, z and Q). When SIDIS datasets are implemented, revisit
# whether SIDIS should also be split into NC and CC, and by hadron type.
SIDIS = _Process(
    "SIDIS",
    "Semi-Inclusive Deep Inelastic Scattering",
    accepted_variables=(_Vars.x, _Vars.z, _Vars.Q2, _Vars.Q, _Vars.y),
    xq2map_function=_dis_xq2map,
    zq2map_function=_sidis_zq2map,
)

SIA = _Process(
    "SIA",
    "Single-Inclusive Annihilation",
    accepted_variables=(_Vars.z, _Vars.xp, _Vars.ph, _Vars.xi, _Vars.sqrts),
    zq2map_function=_sia_zq2map,
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
    "DIJET_3D": DIJET_3D,
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
    "DYP_FT": dataclasses.replace(DY_2L, name="DYP_FT", description="DY Z -> ll (pseudo)rapidity"),
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
    "POS_POL": POS_POL,
    "POS_DY": dataclasses.replace(POS_DIS, name="POS_DY"),
    "INTEG": dataclasses.replace(POS_DIS, name="INTEG", description="Integrability dataset"),
    "PHT": PHT,
    "SINGLETOP": SINGLETOP,
    "SIDIS": SIDIS,
    "SIA": SIA,
}


@Parser
def ValidProcess(process_name) -> Union[_Process, str]:
    return PROCESSES.get(process_name.upper(), process_name.upper())
