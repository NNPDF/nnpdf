"""
This module provides an unique source and definition for all the possible parameters
that a theory card can contain. 
It also implement some utilities function to manage the theory cards. 
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class TheoryCard:
    ID: int # ID number of the theory
    PTO: int # Perturbative order (0 = LO, 1 = NLO, 2 = NNLO ...)
    FNS: str # Flavor number scheme (i.e. FONLL-C)
    DAMP: int # Whether a damping function is applied or not for FONLL
    IC: int # 0 = perturbative charm only , 1 = intrinsic charm allowed
    ModEv: str # DGLAP evolution solution method (EXA or TRN)
    XIR: float # Renormalization scale over the hard scattering scale ratio
    XIF: float # Factorization scale over the hard scattering scale ratio
    NfFF: int # Number of active flavors, only for FFNS or FFN0 schemes
    nfref: int # Number of active flavors at Qref
    Q0: float # [GeV] Parametrization scale
    alphas: float # Value of alpha_s at the scale Qref
    Qref: float # [GeV] Reference scale for alphas and alphaqed
    QED: int # QED correction to strong coupling: 0 = disabled , 1 = allowed
    alphaqed: float # Values of alpha QED at the scale Qref
    HQ: str # Heavy quark mass scheme,  POLE for pole masses (default), MSBAR for running masses (used only in Eko).
    mc: float # [GeV] charm mass
    Qmc: float # [GeV] MSbar mass reference scale of the charm
    kcThr: float # Threshold ratio of the charm
    mb: float # # [GeV] bottom mass
    Qmb: float # [GeV] MSbar mass reference scale of the bottom
    kbThr: float # Threshold ratio of the bottom
    mt: float # # [GeV] top mass
    Qmt: float # [GeV] MSbar mass reference scale of the top
    ktThr: float # Threshold ratio of the top
    CKM: list[float] # CKM matrix elements
    MZ: float # [GeV] Mass of Z
    MW: float # [GeV] Mass of W
    GF: float # Fermi constant
    SIN2TW: float
    TMC: int # Include target mass corrections: 0 = disabled, 1 = leading twist, 2 = higher twist approximated, 3 = higher twist exact
    MP: float # [GeV] Mass of the proton
    Comments: str # Comments on the theory
    IterEv: int = None # Number of iterations for the evolution of the PDF. Defaults to 40 when ModEv = EXA
    ModSV: str = None # Scale variations method in EKO (expanded or exponentiated)
    DAMPPOWERc: int = None # Power of the damping factor in FONLL for the c
    DAMPPOWERb: int = None # Power of the damping factor in FONLL for the b
    n3lo_cf_variation : int = 0 # N3LO coefficient functions variation: -1 = lower bound, 0 = central, 1 = upper bound
    use_fhmruvv: bool = False # N3LO splitting functions approximation: if True use the FHMRUVV parametrization, otherwise use EKO parametrization.
    nf0: int = 3 # Number of active flavors at the parametrization scale Q0