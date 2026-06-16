"""Evolution of n3fit export grids with HOPPET.

This module mirrors the LHAPDF-writing part of :mod:`evolven3fit.evolve`,
but uses HOPPET for the DGLAP evolution instead of applying a precomputed EKO.
"""

from collections import defaultdict
import dataclasses
import functools
import importlib
import logging
import pathlib
import shutil
import sys
import tempfile

from ekobox import genpdf, info_file
import numpy as np

from eko import EKO, basis_rotation
from eko.interpolation import InterpolatorDispatcher
from nnpdf_data import THEORY_CARDS_PATH
from nnpdf_data.theorydbutils import fetch_theory
from validphys.utils import yaml_safe

from . import eko_utils
from .evolve import LOG_FILE, LOGGING_SETTINGS, ExportGrid, dump_evolved_replica

_logger = logging.getLogger(__name__)

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


def _configure_hoppet(hp, theory: HoppetTheory, x_grid, q_values, dy, lnlnq_order, y_order):
    """
    Calling a few hoppet function that need to be called before starting.
    We use hoppet's extended start to make sure we have the same
    range as the lhapdf grid that we use in NNPDF.
    """
    if y_order is not None or lnlnq_order is not None:
        hp.SetYLnlnQInterpOrders(
            -1 if y_order is None else int(y_order), 4 if lnlnq_order is None else int(lnlnq_order)
        )

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


def _evolve_one_exportgrid(
    hp, theory: HoppetTheory, exportgrid: ExportGrid, q_by_nf: list[tuple[float, int]]
) -> dict:
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


def evolve_exportgrid_with_hoppet(
    nnpdf_theory,
    exportgrids: list[ExportGrid],
    dy: float = 0.025,
    y_order: int = None,
    lnlnq_order: int = None,
    eko_path: pathlib.Path = None,
):
    """Same-old, same-old, hoppet version of ``evolve_exportgrid``."""
    # TODO: we start with the same exact check as evolve_exportgrid
    ref = exportgrids[0]
    hessian_fit = ref.hessian
    for egrid in exportgrids:
        assert egrid.q20 == ref.q20, "Different values of q0 found among the exportgrids"
        np.testing.assert_allclose(
            ref.xgrid, egrid.xgrid, err_msg="ExportGrids are not all evaluate at the same x nodes"
        )
        assert (
            hessian_fit == egrid.hessian
        ), "Trying to evolve hessian and non-hessian fit at the same time"
    ######

    # Build the theory and operator cards from the NNPDF theories.
    # _if_ an EKO path is given, extract the operator card directly from there.
    # This ensures that the xgrid/qgrid is exactly the same as the eko
    theory_card, op_card = _build_cards(nnpdf_theory, ref.xgrid, eko_path)
    hoppet_theory = _nnpdf_theory_to_hoppet(nnpdf_theory)
    hoppet_module = _import_hoppet()

    eko_q_by_nf = op_card.raw["mugrid"]
    q_values = sorted({q for q, _ in eko_q_by_nf})
    all_evolved = defaultdict(list)

    try:
        _configure_hoppet(
            hoppet_module, hoppet_theory, ref.xgrid, q_values, dy, lnlnq_order, y_order
        )

        # EKO gets to do the whole set at once by creating one single
        # tensor with replicas as an index.
        # Hoppet instead will do them one by one, but hopefully
        # caches the intermediate results?
        for exportgrid in exportgrids:
            evolved_replica = _evolve_one_exportgrid(
                hoppet_module, hoppet_theory, exportgrid, eko_q_by_nf
            )

            for q, nf in eko_q_by_nf:
                all_evolved[(q**2, nf)].append(evolved_replica[(q, nf)])
    finally:
        hoppet_module.DeleteAll()

    # This is a copy from tehe vanilla evolution but we need the extra np.array
    # to "concatenate" the replicas
    sorted_evolved = {
        key: np.array(value)
        for key, value in sorted(all_evolved.items(), key=lambda item: (item[0][1], item[0][0]))
    }

    info = info_file.build(theory_card, op_card, 1, info_update={})
    info["NumMembers"] = "REPLACE_NREP"
    if hessian_fit:
        info["ErrorType"] = "hessian"
    else:
        info["ErrorType"] = "replicas"
    info["XMin"] = float(ref.xgrid[0])
    info["XMax"] = float(ref.xgrid[-1])
    info["Flavors"] = basis_rotation.flavor_basis_pids
    info.setdefault("NumFlavors", 5)

    # TODO Not sure why this doesn't work ootb, double check
    if info.get("AlphaS_Qs"):
        info["AlphaS_Qs"][0] = info["QMin"]
        info["AlphaS_Qs"][-1] = info["QMax"]

    return info, sorted_evolved


def evolve_exportgrids_into_lhapdf_with_hoppet(
    nnpdf_theory,
    exportgrids: list[ExportGrid],
    output_files: list[pathlib.Path],
    info_path: pathlib.Path,
    finalize: bool = False,
    dy: float = 0.025,
    y_order: int = None,
    lnlnq_order: int = None,
    eko_path: pathlib.Path = None,
):
    """
    Hoppet version of ``evolve_exportgrids_into_lhapdf``
    """
    if len(exportgrids) != len(output_files):
        raise ValueError("The length of output_files and exportgrids must be equal")

    for f in output_files:
        pathlib.Path(f).parent.mkdir(exist_ok=True, parents=True)
    pathlib.Path(info_path).parent.mkdir(exist_ok=True, parents=True)

    # TODO: at this point everything other than _this_ call is a direct copy from the
    # vanilla evolution for development convenience, both should be merged together before merging
    # plz reviewer, scold me if I haven't done so
    info, all_evolved = evolve_exportgrid_with_hoppet(
        nnpdf_theory,
        exportgrids,
        dy=dy,
        y_order=y_order,
        lnlnq_order=lnlnq_order,
        eko_path=eko_path,
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = pathlib.Path(temp_dir)

        if finalize:
            info["NumMembers"] = len(exportgrids)

        genpdf.export.dump_info(temp_path, info)
        temp_info = temp_path / f"{temp_path.stem}.info"
        shutil.move(temp_info, info_path)

        # Dump LHAPDF files as .dat files in blocks of nf
        targetgrid = exportgrids[0].xgrid.tolist()
        q2block_per_nf = defaultdict(list)
        for q2, nf in all_evolved.keys():
            q2block_per_nf[nf].append(q2)

        x_index = {x: idx for idx, x in enumerate(targetgrid)}
        pid_index = {pid: idx for idx, pid in enumerate(info["Flavors"])}

        for enum, (exportgrid, output_file) in enumerate(zip(exportgrids, output_files)):
            replica_idx = exportgrid.replica
            if replica_idx is None and exportgrid.hessian:
                replica_idx = enum
            blocks = []

            for nf, q2grid in q2block_per_nf.items():

                def pdf_xq2(pid, x, Q2):
                    x_idx = x_index[x]
                    pid_idx = pid_index[pid]
                    return x * all_evolved[(Q2, nf)][enum][pid_idx][x_idx]

                block = genpdf.generate_block(
                    pdf_xq2, xgrid=targetgrid, sorted_q2grid=q2grid, pids=info["Flavors"]
                )
                blocks.append(block)

            dat_path = dump_evolved_replica(blocks, temp_path, replica_idx, exportgrid.hessian)
            if not dat_path.exists():
                raise FileNotFoundError(
                    f"The expected {dat_path} file was not found after dumping the blocks"
                )
            shutil.move(dat_path, output_file)


def evolve_fit_with_hoppet(
    fit_folder: pathlib.Path,
    theory_id: int,
    force: bool = False,
    hessian_fit: bool = False,
    dy: float = 0.025,
    y_order: int = None,
    lnlnq_order: int = None,
    eko_path: pathlib.Path = None,
):
    """
    Evolves all the fitted replica in fit_folder/nnfit.
    Follows the same structure as ``evolven3fit.evolve.evolve_fit`` but using HOPPET
    as the backend instead of a pre-computed EKO.
    """
    # Check that hoppet can be imported, otherwise fail early
    _import_hoppet()

    # TODO: at this point this piece is repeated here and in evole_fit
    # so it should be lifted
    fit_folder = pathlib.Path(fit_folder)
    log_path = fit_folder / LOG_FILE
    if log_path.exists():
        if force:
            log_path.unlink()
        else:
            raise FileExistsError(
                f"Log file already exists: {log_file}. Has evolven3fit already been run?"
            )

    log_file = logging.FileHandler(log_path)
    stdout_log = logging.StreamHandler(sys.stdout)
    for log in [log_file, stdout_log]:
        log.setFormatter(LOGGING_SETTINGS["formatter"])
    log_file.setLevel(LOGGING_SETTINGS["level"])
    stdout_log.setLevel(logging.INFO)

    for logger in (_logger, logging.getLogger("evolven3fit")):
        logger.handlers = []
        logger.setLevel(LOGGING_SETTINGS["level"])
        logger.propagate = False
        logger.addHandler(log_file)
        logger.addHandler(stdout_log)
    #############################################

    nnpdf_theory = fetch_theory(THEORY_CARDS_PATH, theory_id, as_dict=False)
    _logger.info("Evolving %s with HOPPET using theory %s", fit_folder, theory_id)

    output_files = []
    exportgrids = []

    for exportgrid_file in fit_folder.glob(f"nnfit/replica_*/{fit_folder.name}.exportgrid"):
        data = yaml_safe.load(exportgrid_file.read_text(encoding="UTF-8"))
        exportgrids.append(ExportGrid(**data, hessian=hessian_fit))
        output_files.append(exportgrid_file.with_suffix(".dat"))

    info_path = fit_folder / "nnfit" / f"{fit_folder.name}.info"
    evolve_exportgrids_into_lhapdf_with_hoppet(
        nnpdf_theory,
        exportgrids,
        output_files,
        info_path,
        dy=dy,
        y_order=y_order,
        lnlnq_order=lnlnq_order,
        eko_path=eko_path,
    )
