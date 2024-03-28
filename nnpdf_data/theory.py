"""
This module provides an unique source and definition for all the possible parameters
that a theory card can contain. 
It also implement some utilities function to manage the theory cards. 
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class _TheoryCard:
    ID: int
    PTO: int
    FNS: str
    DAMP: int
    DAMPPOWERc: None #Power of the damping factor in FONLL for the c
    DAMPPOWERb: None #Power of the damping factor in FONLL for the b
    IC: int
    IB: int
    ModEv: str
    ModSV: None#
    XIR: float
    XIF: float
    NfFF: int
    #MaxNfAs: int
    #MaxNfPdf: int
    nfref: 5
    nf0: 3
    Q0: float
    alphas: float
    Qref: float #TODO: understand if we want to enforce Qref = Qedref and remove Qedref. Note that this is mandatory in eko
    QED: int
    alphaqed: float
    Qedref: float
    SxRes: int#
    SxOrd: str#
    HQ: str
    mc: float
    Qmc: float#
    kcThr: float
    mb: float
    Qmb: float#
    kbThr: float
    mt: float
    Qmt: float#
    ktThr: float
    CKM: list[float]
    MZ: float
    MW: float
    GF: float
    SIN2TW: float#
    TMC: int
    MP: float
    Comments: str
    global_nx: int#
    #EScaleVar: int
    FactScaleVar: None #used in yadism to allow for Fact scale var
    RenScaleVar: None #used in yadism to allow for Ren scale var
    n3lo_cf_variation : None