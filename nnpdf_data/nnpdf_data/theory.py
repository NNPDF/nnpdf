"""
This module provides an unique source and definition for all the possible parameters
that a theory card can contain.
"""

import dataclasses
import logging

DEPRECATED_KEYS = ["MaxNfAs", "SxRes", "SxOrd" "EScaleVar", "Qedref", "global_nx"]

log = logging.getLogger(__name__)


class TheoryCardError(Exception):
    pass


@dataclasses.dataclass(frozen=True)
class TheoryCard:
    ID: int  # ID number of the theory
    PTO: int  # Perturbative order (0 = LO, 1 = NLO, 2 = NNLO ...)
    FNS: str  # Flavor number scheme (i.e. FONLL-C)
    DAMP: int  # Whether a damping function is applied or not for FONLL
    IC: int  # 0 = perturbative charm only , 1 = intrinsic charm allowed
    ModEv: str  # DGLAP evolution solution method (EXA or TRN)
    XIR: float  # Renormalization scale over the hard scattering scale ratio
    XIF: float  # Factorization scale over the hard scattering scale ratio
    NfFF: int  # Number of active flavors, only for FFNS or FFN0 schemes
    QED: int  # QED correction to strong coupling: 0 = disabled , 1 = allowed
    HQ: str = "POLE" # Heavy quark mass scheme,  POLE for pole masses (default), MSBAR for running masses (used currently only in eko).
    mc: float  # [GeV] charm mass
    Qmc: float  # [GeV] MSbar mass reference scale of the charm
    kcThr: float  # Threshold ratio of the charm
    mb: float  # # [GeV] bottom mass
    Qmb: float  # [GeV] MSbar mass reference scale of the bottom
    kbThr: float  # Threshold ratio of the bottom
    mt: float  # # [GeV] top mass
    Qmt: float  # [GeV] MSbar mass reference scale of the top
    ktThr: float  # Threshold ratio of the top
    CKM: list[float]  # CKM matrix elements
    MZ: float  # [GeV] Mass of Z
    MW: float  # [GeV] Mass of W
    GF: float  # Fermi constant
    SIN2TW: float
    TMC: int  # Include target mass corrections: 0 = disabled, 1 = leading twist, 2 = higher twist approximated, 3 = higher twist exact
    MP: float  # [GeV] Mass of the proton
    Comments: str  # Comments on the theory
    MaxNfPdf: int = 5  # Used by pineko and the photon module to define the thresholds
    # Fit theory parameters default
    nf0: int = 4  # Number of active flavors at the parametrization scale Q0
    Q0: float = 1.65  # [GeV] Parametrization scale
    nfref: int = 5  # Number of active flavors at Qref
    Qref: float = 91.2  # [GeV] Reference scale for alphas and alphaqed
    alphas: float = 0.118  # Value of alpha_s at the scale Qref
    alphaqed: float = 0.007496252  # Values of alpha QED at the scale Qref
    # Evolution parameters
    IterEv: int = None  # iterations for the evolution of the PDF. Defaults to 40 when ModEv = EXA
    ModSV: str = None  # Scale variations method in EKO (expanded or exponentiated)
    DAMPPOWERc: int = None  # Power of the damping factor in FONLL for the c
    DAMPPOWERb: int = None  # Power of the damping factor in FONLL for the b
    # N3LO parameters
    n3lo_ad_variation: list = dataclasses.field(
        default_factory=lambda: 7 * [0]
    )  # N3LO anomalous dimension variations
    n3lo_cf_variation: int = (
        0 # N3LO coefficient functions variation: -1 = lower bound, 0 = central, 1 = upper bound
    )
    use_fhmruvv: bool = (
        False # N3LO splitting functions approximation: if True use the FHMRUVV parametrization, otherwise use EKO parametrization.
    )
    ###### Keys for compatibility with old NNPDF theories, their values will be dropped immediately after reading to avoid problems
    ###### they will be set to ``None`` immediately after loading the theory
    MaxNfAs: int = None
    SxRes: int = None
    SxOrd: str = None
    EScaleVar: int = None
    Qedref: float = None
    global_nx: int = None

    def __post_init__(self):
        """Drop deprecated keys and apply some checks"""
        if self.Qedref is not None and self.QED != 0:
            # Check that nobody is trying to use this with a wrong Qedref!
            if self.Qedref != self.Qref:
                self._raise_or_warn(
                    f"Trying to use {self.ID} with {self.Qedref} != {self.Qref}. This is not supported!"
                )

        for key in DEPRECATED_KEYS:
            object.__setattr__(self, key, None)

    def _raise_or_warn(self, msg):
        """Raise an error for new theories and a warning for old ones"""
        if self.ID >= 600:
            raise TheoryCardError(msg)
        log.warning(msg)

    def asdict(self):
        return dataclasses.asdict(self)
