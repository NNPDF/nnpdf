"""
This module provides an unique source and definition for all the possible parameters
that a theory card can contain.
"""

import dataclasses
import logging
from typing import Literal, Optional

DEPRECATED_KEYS = ["MaxNfAs", "SxRes", "SxOrd" "EScaleVar", "Qedref", "global_nx"]

log = logging.getLogger(__name__)


class TheoryCardError(Exception):
    pass


class TheoryNotFoundInDatabase(FileNotFoundError):
    pass


@dataclasses.dataclass(frozen=True)
class TheoryCard:
    ID: int  # ID number of the theory
    PTO: Literal[0, 1, 2, 3]  # Perturbative order (0 = LO, 1 = NLO, 2 = NNLO ...)
    FNS: str  # Flavor number scheme (e.g. FONLL-C)
    IC: int  # 0 = perturbative charm only , 1 = intrinsic charm allowed
    ModEv: Literal["EXA", "TRN"]  # DGLAP evolution solution method (e.g. EXA or TRN)
    XIR: float  # Renormalization scale over the hard scattering scale ratio
    XIF: float  # Factorization scale over the hard scattering scale ratio
    NfFF: int  # Number of active flavors, only for FFNS or FFN0 schemes
    QED: int  # Max order of alpha_qed in the evolution
    Q0: float  # [GeV] Parametrization scale for the fit (and the photon)
    mc: float  # [GeV] charm mass
    Qmc: float  # [GeV] MSbar mass reference scale of the charm
    kcThr: float  # Threshold ratio of the charm
    mb: float  # # [GeV] bottom mass
    Qmb: float  # [GeV] MSbar mass reference scale of the bottom
    kbThr: float  # Threshold ratio of the bottom
    mt: float  # # [GeV] top mass
    Qmt: float  # [GeV] MSbar mass reference scale of the top
    ktThr: float  # Threshold ratio of the top
    # CKM matrix elements (running on the columns first, i.e. CKM_11 is CKM[0] and CKM_12 is CKM[1] and so on)
    CKM: list[float]
    MZ: float  # [GeV] Mass of Z
    MW: float  # [GeV] Mass of W
    GF: float  # Fermi constant
    SIN2TW: float
    TMC: int  # Include target mass corrections: 0 = disabled, 1 = leading twist, 2 = higher twist approximated, 3 = higher twist exact
    MP: float  # [GeV] Mass of the proton
    Comments: str  # Comments on the theory
    # Definition of the reference scale for alpha_s and alpha_qed
    alphas: float  # Value of alpha_s at the scale Qref
    alphaqed: float  # Values of alpha QED at the scale Qref
    Qref: float  # [GeV] Reference scale for alphas and alphaqed
    nfref: Optional[int] = None  # nf at Qref (its default depend on Qref)
    MaxNfPdf: Optional[int] = 5  # Used by pineko and the photon module to define the thresholds
    ## Fit theory parameters default
    # Number of active flavors at the parametrization scale Q0, its default depends on Q0
    nf0: Optional[int] = None
    ## Evolution parameters
    # Heavy quark mass scheme, POLE for pole masses (default), MSBAR for running masses
    HQ: Optional[Literal["POLE", "MSBAR"]] = "POLE"
    # iterations for the evolution of the PDF. Defaults to 60 when ModEv = EXA
    IterEv: Optional[int] = None
    ModSV: Optional[str] = None  # Scale variations method in EKO (expanded or exponentiated)
    ## NFONLL parameters, used only by pineko
    DAMP: Optional[int] = 0  # Whether a damping function is applied or not for FONLL
    DAMPPOWERc: Optional[int] = None  # Power of the damping factor in FONLL for the c
    DAMPPOWERb: Optional[int] = None  # Power of the damping factor in FONLL for the b
    FONLLParts: Optional[str] = None
    PTOEKO: Optional[int] = None
    PTODIS: Optional[int] = None
    ## N3LO parameters
    # N3LO anomalous dimension variations
    n3lo_ad_variation: Optional[list] = dataclasses.field(default_factory=lambda: 7 * [0])
    # N3LO coefficient functions variation: -1 = lower bound, 0 = central, 1 = upper bound
    n3lo_cf_variation: Optional[int] = 0
    # N3LO splitting functions approximation: if True use the FHMRUVV parametrization, otherwise use EKO parametrization.
    use_fhmruvv: Optional[bool] = False
    ###### Keys for compatibility with old NNPDF theories, their values will be dropped immediately after reading to avoid problems
    ###### they will be set to ``None`` immediately after loading the theory
    MaxNfAs: Optional[int] = None
    SxRes: Optional[int] = None
    SxOrd: Optional[str] = None
    EScaleVar: Optional[int] = None
    Qedref: Optional[float] = None
    global_nx: Optional[int] = None

    def __post_init__(self):
        self._set_defaults()
        self._apply_checks()

        for key in DEPRECATED_KEYS:
            object.__setattr__(self, key, None)

    def _apply_checks(self):
        """Apply complex checks to the input keys after parsing is complete"""
        if self.Qedref is not None and self.QED != 0:
            # Check that nobody is trying to use this with a wrong Qedref!
            if self.Qedref != self.Qref:
                self._raise_or_warn(
                    f"Trying to use {self.ID} with {self.Qedref} != {self.Qref}. This is not supported!"
                )

        if self.XIF == 1 and self.ModSV is not None:
            raise TheoryCardError(
                f"Theory: {self.ID}, error: XIF is {self.XIF} while ModSV is {self.ModSV}. "
                "If XIF is equal to 1.0, ModSV should not be defined."
            )
        elif self.XIF != 1 and self.ModSV not in ["expanded", "exponentiated"]:
            raise TheoryCardError(
                f"Theory: {self.ID}, error: XIF is {self.XIF} while ModSV is {self.ModSV}. "
                "If XIF is different from 1.0, ModSV should be either 'expanded' or 'exponentiated'."
            )

        if self.DAMP != 0 and "FONLL" in self.FNS:
            # Check the damp powers are being used
            if self.DAMPPOWERb is None:
                self._raise_or_warn("DAMPOWERb key needs to be given when DAMP and FONLL are used")
            if self.DAMPPOWERc is None and self.IC == 0:
                self._raise_or_warn("DAMPOWERc key needs to be given when DAMP and FONLL are used")

    def find_nf(self, mu):
        """Given a value of q, find the corresponding number of flavours
        for this theory."""
        if mu < self.mc * self.kcThr:
            return 3
        elif mu < self.mb * self.kbThr:
            return 4
        elif mu < self.mt * self.ktThr:
            return 5
        return 6

    def _set_defaults(self):
        """Function to be called by __post_init__ to set defaults that might depends
        on other variables"""
        if self.ModEv == "EXA" and self.IterEv is None:
            object.__setattr__(self, "IterEv", 60)
        if self.ModEv == "TRN" and self.IterEv is not None:
            raise TheoryCardError(
                f"IterEv should not be given when ModEv=TRN, received {self.IterEv}"
            )
        if self.nf0 is None:
            object.__setattr__(self, "nf0", self.find_nf(self.Q0))
        if self.nfref is None:
            object.__setattr__(self, "nfref", self.find_nf(self.Qref))

    def _raise_or_warn(self, msg):
        """Raise an error for new theories and a warning for old ones"""
        if self.ID >= 600:
            raise TheoryCardError(f"Theory: {self.ID}, error: {msg}")
        log.warning(msg)

    def asdict(self):
        return dataclasses.asdict(self)
