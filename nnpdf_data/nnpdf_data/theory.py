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
    DAMPPOWERc: int = None # Power of the damping factor in FONLL for the c
    DAMPPOWERb: int = None # Power of the damping factor in FONLL for the b
    IC: int # 0 = perturbative charm only , 1 = intrinsic charm allowed
    ModEv: str # DGLAP evolution solution method (EXA or TRN)
    IterEv: int = None # Number of iterations for the evolution of the PDF. Defaults to 40 when ModEv = EXA
    ModSV: str = None # Scale variations method in EKO (expanded or exponentiated)
    XIR: float # Renormalization scale ratio
    XIF: float # Factorization scale ratio
    NfFF: int # Number of active flavors
    nfref: 5 # Number of active flavors at Qref
    nf0: 3 # Number of active flavors at the parametrization scale Q0
    Q0: float # Parametrization scale
    alphas: float # Value of alpha_s at the scale Qref
    Qref: float # Reference scale for alphas and alphaqed
    QED: int # Whether QED effects are taken into account
    alphaqed: float # Values of alpha QED at the scale Qref
    HQ: str # Heavy quark mass scheme,  POLE for pole masses (default), MSBAR for running masses (used only in Eko).
    mc: float # Pole mass of the charm
    Qmc: float # MSbar mass reference scale of the charm
    kcThr: float # Threshold ratio of the charm
    mb: float # Pole mass of the bottom
    Qmb: float # MSbar mass reference scale of the bottom
    kbThr: float # Threshold ratio of the bottom
    mt: float # Pole mass of the top
    Qmt: float # MSbar mass reference scale of the top
    ktThr: float # Threshold ratio of the top
    CKM: list[float] # CKM matrix elements
    MZ: float # Mass of Z
    MW: float # Mass of W
    GF: float # Fermi constant
    SIN2TW: float
    TMC: int # Time like
    MP: float # Mass of the proton
    Comments: str # Comments on the theory
    FactScaleVar: None #used in yadism to allow for Fact scale var
    RenScaleVar: None #used in yadism to allow for Ren scale var
    n3lo_cf_variation : None # Variation of the N3LO coefficient function approximation
    use_fhmruvv: None