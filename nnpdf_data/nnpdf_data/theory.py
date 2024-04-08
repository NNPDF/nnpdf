"""
This module provides an unique source and definition for all the possible parameters
that a theory card can contain. 
It also implement some utilities function to manage the theory cards. 
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class TheoryCard:
    ID: int
    PTO: int
    FNS: str
    DAMP: int
    DAMPPOWERc: int = None #Power of the damping factor in FONLL for the c
    DAMPPOWERb: int = None #Power of the damping factor in FONLL for the b
    IC: int
    ModEv: str # EXA or TRN
    IterEv: int = None # Number of iterations for the evolution of the PDF. Defaults to 40 when ModEv = EXA
    ModSV: str = None
    XIR: float
    XIF: float
    NfFF: int
    nfref: 5
    nf0: 3
    Q0: float
    alphas: float
    Qref: float
    QED: int
    alphaqed: float
    HQ: str
    mc: float
    Qmc: float
    kcThr: float
    mb: float
    Qmb: float
    kbThr: float
    mt: float
    Qmt: float
    ktThr: float
    CKM: list[float]
    MZ: float
    MW: float
    GF: float
    SIN2TW: float
    TMC: int
    MP: float
    Comments: str
    FactScaleVar: None #used in yadism to allow for Fact scale var
    RenScaleVar: None #used in yadism to allow for Ren scale var
    n3lo_cf_variation : None
    use_fhmruvv: None