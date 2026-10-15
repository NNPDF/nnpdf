#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
postfit - fit output processing tool

postfit constructs a final PDF fit ensemble and LHAPDF set from the
output of a batch of `nnfit` jobs. The jobs are assessed on the basis
of standard NNPDF fit veto criteria as described in the validphys
fitveto module. Passing replicas are symlinked to the [results]/postfit
directory in the standard NNPDF layout and also in the standard LHAPDF
layout.
"""
__authors__ = 'Nathan Hartland, Zahari Kassabov'

import argparse
from glob import glob
import itertools
import json
import logging
import os.path
import pathlib
import re
import shutil
import sys

import numpy as np

from reportengine import colors
from validphys import fitdata, fitveto, lhio, replica_checks
from validphys.core import PDF
from validphys.fitveto import INTEG_THRESHOLD, NSIGMA_DISCARD_ARCLENGTH, NSIGMA_DISCARD_CHI2

# With per-replica checks, replicas failing a veto can complete the set only if their chi2 is at
# most this factor times the mean chi2 of the replicas passing all vetoes
FILL_CHI2_FACTOR = 2.0
from validphys.lhaindex import paths_prepend
from validphys.lhapdf_compatibility import make_pdf
from validphys.loader import Loader
from validphys.utils import tempfile_cleaner, yaml_safe

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(colors.ColorHandler())


def splash():
    print("                                                        ")
    print("  ██████╗  ██████╗ ███████╗████████╗███████╗██╗████████╗")
    print("  ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝██╔════╝██║╚══██╔══╝")
    print("  ██████╔╝██║   ██║███████╗   ██║   █████╗  ██║   ██║   ")
    print("  ██╔═══╝ ██║   ██║╚════██║   ██║   ██╔══╝  ██║   ██║   ")
    print("  ██║     ╚██████╔╝███████║   ██║   ██║     ██║   ██║   ")
    print("  ╚═╝      ╚═════╝ ╚══════╝   ╚═╝   ╚═╝     ╚═╝   ╚═╝   ")
    print("  __coredevs__: S.C, N.H. Z.K.\n")


def relative_symlink(source, dest):
    """Forms a relative symlink between 'source' and 'dest'"""
    relativepath = os.path.relpath(source, dest.parents[0])
    os.symlink(relativepath, dest)


def set_lhapdf_info(info_path, nrep):
    """Sets the LHAPDF info file NumMembers field"""
    with open(info_path, 'r+') as f:
        txt = f.read()
        f.seek(0)
        f.write(txt.replace('REPLACE_NREP', str(nrep + 1)))
        f.truncate()


class PostfitError(Exception):
    """Exception raised when postfit cannot succeed and knows why"""

    pass


class FatalPostfitError(Exception):
    """Exception raised when some corrupted input is detected"""

    pass


def per_replica_positivity(valid_paths, fitname, runcard, postfit_path):
    """Per-replica positivity mask (see :py:mod:`validphys.replica_checks`).

    Replicas whose ``.json`` was written with a per-replica ``pos_state`` (n3fit run with
    ``parameters::bayesian::per_replica::positivity: true``) use it directly; for the others
    the positivity is recomputed from their ``.exportgrid``."""
    mask = np.zeros(len(valid_paths), dtype=bool)
    to_compute = []
    for i, path in enumerate(valid_paths):
        info = json.loads((pathlib.Path(path) / f"{fitname}.json").read_text(encoding="utf-8"))
        if info.get(replica_checks.PER_REPLICA_JSON_KEY, False):
            mask[i] = info["pos_state"] == "POS_PASS"
        else:
            to_compute.append(i)
    log.info(
        f"Per-replica positivity: {len(valid_paths) - len(to_compute)} replicas checked by n3fit, "
        f"{len(to_compute)} recomputed from their exportgrid"
    )
    if to_compute:
        passing, losses = replica_checks.replicas_positivity(
            [valid_paths[i] for i in to_compute], fitname, runcard
        )
        mask[to_compute] = passing
        worst = losses.max(axis=0).sort_values(ascending=False)
        log.info("Largest positivity loss per set over recomputed replicas:\n%s", worst.to_string())
        losses.index = [pathlib.Path(valid_paths[i]).name for i in to_compute]
        losses.insert(0, "passes", passing)
        losses.to_csv(postfit_path / "positivity_losses.csv")
    log_positivity_counts(valid_paths, mask, runcard)
    return mask


def bnn_samples_per_model(runcard):
    """n_bnn_samples of a BNN fit runcard (``parameters::bayesian`` or directly under
    ``parameters``), None otherwise"""
    parameters = (runcard or {}).get("parameters", {})
    return (parameters.get("bayesian") or {}).get("n_bnn_samples", parameters.get("n_bnn_samples"))


def spread_over_models(paths, nrep, runcard):
    """The first ``nrep`` of ``paths``; for a BNN fit with several trained models
    (replicas ``[m*n_bnn_samples, (m+1)*n_bnn_samples)`` of model m), the replicas are taken
    in turn from every model (in replica order within a model), so that every model is
    represented instead of only the first ones."""
    n_samples = bnn_samples_per_model(runcard)
    if not n_samples:
        return paths[:nrep]
    by_model = {}
    for path in paths:
        by_model.setdefault(int(pathlib.Path(path).name.split("_")[-1]) // n_samples, []).append(path)
    if len(by_model) == 1:
        return paths[:nrep]
    selected = [p for group in itertools.zip_longest(*by_model.values()) for p in group if p]
    log.info(
        f"Selecting {nrep} of {len(paths)} passing replicas in turn from the {len(by_model)} "
        f"trained models"
    )
    return selected[:nrep]


def log_positivity_counts(valid_paths, mask, runcard):
    """Log how many replicas fully respect positivity and, for a BNN fit with several trained
    models (n_bnn_samples replicas each), how the passes are spread over the models."""
    n_pass = int(mask.sum())
    log.info(
        f"Per-replica positivity: {n_pass} of {len(mask)} replicas fully respect positivity, "
        f"{len(mask) - n_pass} violate it"
    )
    n_samples = bnn_samples_per_model(runcard)
    if not n_samples or n_samples >= len(mask):
        return
    model = np.array([int(pathlib.Path(path).name.split("_")[-1]) // n_samples for path in valid_paths])
    passes = np.array([mask[model == m].mean() for m in np.unique(model)])
    log.info(
        f"Per-replica positivity per trained model ({len(passes)} models x {n_samples} samples): "
        f"pass fraction min {passes.min():.2f}, median {np.median(passes):.2f}, max {passes.max():.2f}; "
        f"{int(np.sum(passes == 0))} models with no passing sample, "
        f"{int(np.sum(passes == 1))} models with all samples passing"
    )


def filter_replicas(
    postfit_path,
    nnfit_path,
    fitname,
    chi2_threshold,
    arclength_threshold,
    integ_threshold,
    per_replica_checks=False,
    runcard=None,
):
    """Find the paths of all replicas passing the standard NNPDF fit vetoes
    as defined in fitveto.py. Returns a list of the replica directories that pass.
    If ``per_replica_checks`` is True, a per-replica positivity veto is added."""
    # This glob defines what is considered a valid replica
    # all the following code uses paths from this glob
    # We sort the paths so that the selection of replicas is deterministic
    all_replicas = sorted(glob(f"{nnfit_path}/replica_*/"))
    valid_paths = [path for path in all_replicas if fitdata.check_replica_files(path, fitname)]
    log.info(f"{len(all_replicas)} total replicas found")
    log.info(f"{len(valid_paths)} valid replicas found")

    if len(valid_paths) == 0:
        raise PostfitError("No valid replicas found")

    # Read FitInfo and compute vetoes
    fitinfo = []
    for path in valid_paths:
        try:
            fitinfo.append(fitdata.load_fitinfo(pathlib.Path(path), fitname))
        except Exception as e:
            raise FatalPostfitError(
                f"Corrupted replica replica at {path}. "
                f"Error when loading replica information:\n {e}"
            ) from e
    positivity_mask = None
    if per_replica_checks:
        positivity_mask = per_replica_positivity(valid_paths, fitname, runcard, postfit_path)
    fit_vetoes = fitveto.determine_vetoes(
        fitinfo, chi2_threshold, arclength_threshold, integ_threshold, positivity_mask
    )
    fitveto.save_vetoes_info(
        fit_vetoes,
        chi2_threshold,
        arclength_threshold,
        integ_threshold,
        postfit_path / "veto_count.json",
    )

    for key in fit_vetoes:
        log.info("%d replicas pass %s" % (sum(fit_vetoes[key]), key))
    passing_paths = list(itertools.compress(valid_paths, fit_vetoes["Total"]))

    # With per-replica checks: the other valid replicas with chi2 <= FILL_CHI2_FACTOR times the
    # mean chi2 of the passing ones, to complete the set if too few pass
    fill_paths = []
    if positivity_mask is not None:
        chi2 = np.array([info.chi2 for info in fitinfo])
        max_chi2 = FILL_CHI2_FACTOR * chi2[fit_vetoes["Total"]].mean()
        eligible = ~fit_vetoes["Total"] & (chi2 <= max_chi2)
        fill_paths = list(itertools.compress(valid_paths, eligible))
        log.info(
            f"{len(fill_paths)} of the {int(np.sum(~fit_vetoes['Total']))} other valid replicas have "
            f"chi2 <= {FILL_CHI2_FACTOR:g} x the passing mean ({max_chi2:.3f}) and can complete the set"
        )
    return passing_paths, fill_paths


def type_fitname(fitname: str):
    """Ensure the sanity of the fitname"""
    fitpath = pathlib.Path(fitname).absolute()
    # Accept only [a-Z, 0-9, -, _]
    sane_name = re.compile(r'[\w\-]+')
    if sane_name.fullmatch(fitpath.name) is None:
        new_name = "-".join(i for i in sane_name.findall(fitname))
        raise argparse.ArgumentTypeError(
            "Only alphanumeric characters, _ and - are allowed in the fit name. "
            f"Please, re-run postfit after renaming the fit: ~$ vp-fitrename {fitpath} {new_name}"
        )
    return fitpath


def _postfit(
    results: str,
    nrep: int,
    chi2_threshold: float,
    arclength_threshold: float,
    integ_threshold: float,
    per_replica_checks: bool = False,
):
    result_path = pathlib.Path(results).resolve()
    fitname = result_path.name

    # Per-replica checks: on if requested on the command line or in the fit runcard
    runcard_path = result_path / "filter.yml"
    runcard = yaml_safe.load(runcard_path.read_text(encoding="utf-8")) if runcard_path.exists() else {}
    if replica_checks.postfit_veto_requested(runcard):
        per_replica_checks = True
    if per_replica_checks and not runcard:
        raise PostfitError(f"Per-replica checks need the fit runcard, not found at {runcard_path}")

    # Paths
    nnfit_path = result_path / 'nnfit'  # Path of nnfit replica output
    final_postfit_path = result_path / 'postfit'
    # Create a temporary path to store work in progress and move it to
    # the final location in the end,
    with tempfile_cleaner(
        root=result_path,
        exit_func=shutil.move,
        exc=(KeyboardInterrupt, PostfitError),
        prefix="postfit_work_deleteme_",
        dst=final_postfit_path,
    ) as postfit_path:

        LHAPDF_path = postfit_path / fitname  # Path for LHAPDF grid output

        if not fitdata.check_nnfit_results_path(result_path):
            raise PostfitError('Postfit cannot find a valid results path')
        if not fitdata.check_lhapdf_info(result_path, fitname):
            raise PostfitError(
                'Postfit cannot find a valid LHAPDF info file\n'
                'Note: Perhaps you forgot to run evolven3fit? See\n'
                'https://docs.nnpdf.science/tutorials/run-fit.html#running-the-fitting-code'
            )

        nrep = int(nrep)
        log.warning(f"Postfit aiming for {nrep} replicas")

        # Generate postfit and LHAPDF directory
        if final_postfit_path.is_dir():
            log.warning(f"Removing existing postfit directory: {final_postfit_path}")
            shutil.rmtree(final_postfit_path)
        os.mkdir(LHAPDF_path)

        # Setup postfit log
        postfitlog = logging.FileHandler(postfit_path / 'postfit.log', mode='w')
        log.addHandler(postfitlog)

        # Perform postfit selection
        if per_replica_checks:
            log.warning("Per-replica positivity checks are on")
        passing_paths, fill_paths = filter_replicas(
            postfit_path,
            nnfit_path,
            fitname,
            chi2_threshold,
            arclength_threshold,
            integ_threshold,
            per_replica_checks=per_replica_checks,
            runcard=runcard,
        )
        if len(passing_paths) >= nrep:
            # Select nrep passing replicas: the first ones, or for a BNN fit with several
            # trained models taken in turn from every model
            selected_paths = spread_over_models(passing_paths, nrep, runcard)
        else:
            # With per-replica checks: complete the set with other valid replicas of acceptable
            # chi2; if there are not enough, continue with all of them
            n_missing = nrep - len(passing_paths)
            if per_replica_checks and len(fill_paths) < n_missing:
                log.warning(
                    f"Postfit requested {nrep} replicas, but only {len(passing_paths)} pass all vetoes "
                    f"and {len(fill_paths)} more have chi2 <= {FILL_CHI2_FACTOR:g} x the passing mean: "
                    f"continuing with {len(passing_paths) + len(fill_paths)} replicas"
                )
                n_missing = len(fill_paths)
            if len(fill_paths) < n_missing:
                raise PostfitError(
                    f"Number of requested replicas is too large: {len(passing_paths)} pass all "
                    f"vetoes and {len(fill_paths)} other valid replicas"
                )
            added = fill_paths[:n_missing]
            log.warning(
                f"Only {len(passing_paths)} replicas pass all vetoes: adding {n_missing} other "
                f"valid replicas (listed in replicas_added.txt)"
            )
            (postfit_path / "replicas_added.txt").write_text(
                "\n".join(pathlib.Path(path).name for path in added) + "\n", encoding="utf-8"
            )
            selected_paths = passing_paths + added

        # Copy info file
        info_source_path = nnfit_path.joinpath(f'{fitname}.info')
        info_target_path = LHAPDF_path.joinpath(f'{fitname}.info')
        shutil.copy2(info_source_path, info_target_path)
        set_lhapdf_info(info_target_path, len(selected_paths))

        # Generate symlinks
        for drep, source_path in enumerate(selected_paths, 1):
            # Symlink results to postfit directory
            source_dir = pathlib.Path(source_path).resolve()
            target_dir = postfit_path.joinpath('replica_%d' % drep)
            relative_symlink(source_dir, target_dir)
            # Symlink results to pdfset directory
            source_grid = source_dir.joinpath(fitname + '.dat')
            target_file = f'{fitname}_{drep:04d}.dat'
            target_grid = LHAPDF_path.joinpath(target_file)
            relative_symlink(source_grid, target_grid)

        log.info(f"{len(selected_paths)} replicas written to the postfit folder")

        # Generate final PDF with replica 0
        log.info("Beginning construction of replica 0")
        # It's important that this is prepended, so that any existing instance of
        # `fitname` is not read from some other path
        paths_prepend(postfit_path)
        generatingPDF = PDF(fitname)
        lhio.generate_replica0(generatingPDF)

        # Test replica 0
        try:
            make_pdf(fitname, 0)
        except RuntimeError as e:
            raise PostfitError("CRITICAL ERROR: Failure in reading replica zero") from e

    log.info("\n\n*****************************************************************\n")
    log.info("Postfit complete")
    log.info(
        f"""
    If you have access to the NNPDF server, please upload your results with
        ~$ vp-upload {result_path}
    and install with
        ~$ vp-get fit {fitname}

    Otherwise, you can manually install the fit as a validphys result by copying or linking
        {result_path}
    into
        {Loader().resultspath}
    e.g.,
        ~$ ln -s {result_path} {Loader().resultspath}
"""
    )
    log.info("*****************************************************************\n\n")


def main():
    splash()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('nrep', type=int, help="Number of desired replicas")
    parser.add_argument(
        "result_path", type=type_fitname, help="Folder containing the results of the fit"
    )
    parser.add_argument(
        '--chi2-threshold',
        nargs='?',
        default=NSIGMA_DISCARD_CHI2,
        help="The number of standard deviations in the chi2, calculated over PDF replicas, "
        f" above which the replicas are cut. The default is {NSIGMA_DISCARD_CHI2}.",
        type=float,
    )
    parser.add_argument(
        '--arclength-threshold',
        nargs='?',
        default=NSIGMA_DISCARD_ARCLENGTH,
        help="The number of standard deviations in the arclength, calculated over PDF replicas, "
        f" above which the replicas are cut. The default is {NSIGMA_DISCARD_ARCLENGTH}.",
        type=float,
    )
    parser.add_argument(
        '--integrability-threshold',
        nargs='?',
        default=INTEG_THRESHOLD,
        help="The maximum value allowed for integrable distributions at small-x. "
        f"The default is {INTEG_THRESHOLD}.",
        type=float,
    )
    parser.add_argument(
        '--per-replica-checks',
        action='store_true',
        help="Veto replicas failing positivity individually (see validphys.replica_checks). "
        "Switched on automatically if the fit runcard sets "
        "parameters::bayesian::per_replica::postfit_veto. "
        "Use this flag for fits that were run without it. If fewer than nrep replicas pass, "
        f"the set is completed with other valid replicas of chi2 <= {FILL_CHI2_FACTOR:g} x the "
        "passing mean; if there are not enough, postfit continues with all of those.",
    )
    parser.add_argument('-d', '--debug', action='store_true', help='show debug messages')
    args = parser.parse_args()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    try:
        _postfit(
            args.result_path,
            args.nrep,
            args.chi2_threshold,
            args.arclength_threshold,
            args.integrability_threshold,
            args.per_replica_checks,
        )
    except PostfitError as e:
        log.error(f"Error in postfit:\n{e}")
        sys.exit(1)
    except FatalPostfitError as e:
        log.error(f"Corrupted input encountered")
        print(colors.color_exception(e.__class__, e, e.__traceback__), file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        log.critical(f"Bug in postfit occurred. Please report it.")
        print(colors.color_exception(e.__class__, e, e.__traceback__), file=sys.stderr)
        sys.exit(1)
