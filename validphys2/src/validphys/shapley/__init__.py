"""Shapley value analysis for NNPDF PDF flavour importance.

Subpackage of validphys that provides tools for computing exact Shapley
values to quantify the importance of individual PDF flavours in
constraining experimental observables.

Uses the external ``shapley_values`` package for the problem-agnostic
Shapley computation, and wraps NNPDF-specific components:
  - Dense FK tensor convolution for fast repeated chi2 evaluation
  - Gaussian perturbation of PDF grid values
  - MSR + VSR sum rule enforcement
  - Evolution / physical-flavor basis support

Typical usage via CLI::

    vp-shapley runcards/sv_dis_hera.yaml

Or from Python::

    from validphys.shapley import (
        setup_observables,
        NNPDFShapleyAnalyzer,
    )
"""

from .setup import (
    setup_observables,
    get_pdf_grid_values,
    get_pdf_flavor_grid_values,
    get_pdf_grid_values_all14,
    FKEntry,
    NNPDFObservable,
    FLAVOR_PDG_NAMES,
)
from .perturbation import (
    gaussian_profile,
    apply_gaussian_perturbation,
    PERTURBATION_MODES,
    PERTURBATION_XSPACES,
)
from .sumrules import (
    gen_integration_input,
    compute_sumrule_normalization,
)
from .analyzer import NNPDFShapleyAnalyzer

__all__ = [
    # Setup
    "setup_observables",
    "get_pdf_grid_values",
    "get_pdf_flavor_grid_values",
    "get_pdf_grid_values_all14",
    "FKEntry",
    "NNPDFObservable",
    "FLAVOR_PDG_NAMES",
    # Perturbation
    "gaussian_profile",
    "apply_gaussian_perturbation",
    "PERTURBATION_MODES",
    "PERTURBATION_XSPACES",
    # Sum rules
    "gen_integration_input",
    "compute_sumrule_normalization",
    # Analyzer
    "NNPDFShapleyAnalyzer",
]
