"""Evolution of n3fit export grids with HOPPET.

This module mirrors the LHAPDF-writing part of :mod:`evolven3fit.evolve`,
but uses HOPPET for the DGLAP evolution instead of applying a precomputed EKO.
"""

from collections import defaultdict
import dataclasses
import functools
import importlib
import pathlib
import tempfile

from ekobox import info_file
import numpy as np

from eko import EKO, basis_rotation
from eko.interpolation import InterpolatorDispatcher
from nnpdf_data import THEORY_CARDS_PATH
from nnpdf_data.theorydbutils import fetch_theory

from . import eko_utils

HOPPET_QCD_PIDS = (-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6)


def _import_hoppet():
    """Import hoppet lazily since it is not part of the default installation."""
    try:
        return importlib.import_module("hoppet")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError("Can't evolve with hoppet without having hoppet installed") from e


@dataclasses.dataclass(frozen=True)
class HoppetTheory:
    """HOPPET settings translated from an NNPDF theory card."""

    alphas: float
    qref: float
    q0: float
    nloop: int
    mur_over_q: float
    masses: tuple[float, float, float]
    mass_scheme: str
    fns: str
    max_nf_pdf: int


# TODO:
# Hoppet asks for a function with a signature (x, Q), in this case Q would be fixed
# so we need to create an interpolator in x following the exportgrid.
# It would be nice if
#   a) The exportgrid itself offered the interpolation (this is useful!)
#   b) We had some very fine interpolation ready to go
#   c) Pass directly the output of the NN!!
#   d) This could also be done by eko
#   e) So, on second thoughts, we want an exportgrid-compatible NN-interface :) # TODO
# we could
def _make_hoppet_pdf_callback(exportgrid, interpolation_degree=4):
    """Creates an interpolation pdf(x) given an exportgrid.
    NOTE: the function includes Q but this is dropped anyway
    """
    xgrid = exportgrid.xgrid
    pid_columns = [exportgrid.pids.index(pid) for pid in HOPPET_QCD_PIDS]
    pdfgrid = exportgrid.pdfgrid[:, pid_columns]

    dispatcher = InterpolatorDispatcher(xgrid, interpolation_degree, mode_N=False)

    @functools.cache  # assuming that an infinite cache is safe here
    def pdf_callback(x, _q):
        x = float(x)
        if x <= xgrid[0]:
            return pdfgrid[0].tolist()
        if x >= xgrid[-1]:
            return pdfgrid[-1].tolist()
        weights = dispatcher.get_interpolation([x])[0]
        return (weights @ pdfgrid).tolist()

    return pdf_callback


def _nnpdf_theory_to_hoppet(nnpdf_theory, polarized=False):
    """Translation layer between the NNPDF theory card and hoppet's parameters."""
    if polarized:
        raise NotImplementedError("Polarized not implemented for hoppet")

    if nnpdf_theory.QED != 0:
        raise NotImplementedError("QED evolution not yet supported for hoppet")

    if nnpdf_theory.PTO > 2:
        raise NotImplementedError("Only up to NNLO for now with hoppet")

    masses = (nnpdf_theory.mc, nnpdf_theory.mb, nnpdf_theory.mt)
    # Note: for hoppet this is both the masses and the thresholds

    nloop = nnpdf_theory.PTO + 1
    fns = nnpdf_theory.FNS.split("-")[0]

    return HoppetTheory(
        alphas=nnpdf_theory.alphas,
        qref=nnpdf_theory.Qref,
        q0=nnpdf_theory.Q0,
        nloop=nloop,
        mur_over_q=nnpdf_theory.XIR,
        masses=masses,
        mass_scheme=nnpdf_theory.HQ,
        fns=fns,
        max_nf_pdf=nnpdf_theory.MaxNfPdf,
    )


def _build_cards(nnpdf_theory, x_grid, eko_path: pathlib.Path = None):
    """Build the operator and theory cards that will later be used to define the
    operators by hoppet. If an EKO is given, use that for the operator card instead."""
    theory_card, op_card = eko_utils.construct_eko_cards(nnpdf_theory.asdict(), x_grid)
    if eko_path is None:
        return theory_card, op_card

    # If we have an eko, load it and read the operator card
    # this is a EKO-sized GBs penalty on /tmp, isn't there a get-eko-metadata function??
    # we should clean after oulseves quickly
    with tempfile.TemporaryDirectory() as temp_dir:
        op_card_eko = EKO.read(eko_path, dest=pathlib.Path(temp_dir)).operator_card
    return theory_card, op_card_eko


def _configure_hoppet(hp, theory: HoppetTheory, x_grid, q_values, dy=0.025):
    """
    Calling a few hoppet function that need to be called before starting.
    We use hoppet's extended start to make sure we have the same
    range as the lhapdf grid that we use in NNPDF.
    """
    hp.SetExactDGLAP(True, True)

    xmin = float(x_grid[0])
    ymax = -np.log(xmin)

    qmin = min(np.min(q_values), theory.q0)
    qmax = max(np.max(q_values), theory.q0)
    hp.StartExtended(
        ymax,
        dy,
        qmin,
        qmax,
        dy / 4.0,  # hoppet's suggestion
        theory.nloop,
        -6,  # interpolation order
        hp.factscheme_MSbar,
    )

    masses = list(theory.masses)
    if theory.max_nf_pdf < 5:
        masses[1] = max(qmax * 2.0, masses[1])
    if theory.max_nf_pdf < 6:
        masses[2] = max(qmax * 2.0, masses[2])

    hp.SetPoleMassVFN(*masses)
    if theory.fns == "FFNS":
        hp.SetFFN(theory.max_nf_pdf)
    elif theory.mass_scheme == "MSBAR":
        hp.SetMSbarMassVFN(*masses)
    else:
        hp.SetPoleMassVFN(*masses)


def _nudge_q_at_threshold(q, nf, q_by_nf):
    """Nudge Q below/above threshold for the changes of NF."""
    matching_nfs = [other_nf for other_q, other_nf in q_by_nf if np.isclose(other_q, q)]
    if len(matching_nfs) < 2:
        return q

    if nf == min(matching_nfs):
        return np.nextafter(q, q - 1.0)
    elif nf == max(matching_nfs):
        return np.nextafter(q, q + 1.0)
    else:
        return NotImplementedError("If you are playing with a 3-way threshold you can modify this")


def _evolve_one_exportgrid(hp, theory, exportgrid, q_by_nf):
    """
    Receives one single exportgrid, and calls hoppet-evolve on it.
    Note, at this point hoppet must have already been started.
    """
    hp.Evolve(
        theory.alphas,
        theory.qref,
        theory.nloop,
        theory.mur_over_q,
        _make_hoppet_pdf_callback(exportgrid),
        np.sqrt(exportgrid.q20),
    )

    evolved = {}
    for q, nf in q_by_nf:
        q_eval = _nudge_q_at_threshold(q, nf, q_by_nf)
        evolved[(q, nf)] = _eval_pdf_grid(hp, exportgrid.xgrid, q_eval)
    return evolved


def _eval_pdf_grid(hp, x_grid, q):
    """Evaluate the PDF grid at the given value of Q.
    This output is (flavour, x) to match what the rest of the evolution functions expect.
    """
    values = []
    for x in x_grid:
        xpdf = hp.Eval(x, q)
        # This can be avoided by having the translation layer ready, later
        xpdf_dict = dict(zip(HOPPET_QCD_PIDS, xpdf))
        xpdf_dict[22] = 0.0
        values.append([xpdf_dict[i] / x for i in basis_rotation.flavor_basis_pids])

    return np.array(values).T


def evolve_exportgrids_with_hoppet(eko_path, exportgrids, theory_id):
    """ """
    ref = exportgrids[0]
    nnpdf_theory = fetch_theory(THEORY_CARDS_PATH, theory_id, as_dict=False)
    theory_card, op_card = _build_cards(nnpdf_theory, ref.xgrid, eko_path)

    # Build the theory and operator cards from the NNPDF theories.
    # _if_ an EKO path is given, extract the operator card directly from there.
    # This ensures that the xgrid/qgrid is exactly the same as the eko
    hoppet_theory = _nnpdf_theory_to_hoppet(nnpdf_theory)
    hoppet_module = _import_hoppet()

    eko_q_by_nf = op_card.raw["mugrid"]
    q_values = sorted({q for q, _ in eko_q_by_nf})
    all_evolved = defaultdict(list)

    try:
        _configure_hoppet(hoppet_module, hoppet_theory, ref.xgrid, q_values)

        for exportgrid in exportgrids:
            evolved_replica = _evolve_one_exportgrid(
                hoppet_module, hoppet_theory, exportgrid, eko_q_by_nf
            )

            for q, nf in eko_q_by_nf:
                all_evolved[(q**2, nf)].append(evolved_replica[(q, nf)])
    finally:
        hoppet_module.DeleteAll()

    info = info_file.build(theory_card, op_card, 1, info_update={})
    sorted_evolved = {
        key: np.array(value)
        for key, value in sorted(all_evolved.items(), key=lambda item: (item[0][1], item[0][0]))
    }
    return info, sorted_evolved
