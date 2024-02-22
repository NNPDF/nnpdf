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
import logging
import os.path
import pathlib
import re
import shutil
import sys

import lhapdf

from reportengine import colors
from validphys import fitdata, fitveto, lhio
from validphys.core import PDF
from validphys.fitveto import INTEG_THRESHOLD, NSIGMA_DISCARD_ARCLENGTH, NSIGMA_DISCARD_CHI2
from validphys.utils import tempfile_cleaner

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


def filter_replicas(
    postfit_path, nnfit_path, fitname, chi2_threshold, arclength_threshold, integ_threshold
):
    """Find the paths of all replicas passing the standard NNPDF fit vetoes
    as defined in fitveto.py. Returns a list of the replica directories that pass."""
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
    fit_vetoes = fitveto.determine_vetoes(
        fitinfo, chi2_threshold, arclength_threshold, integ_threshold
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
    return passing_paths


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
    at_least_nrep: bool,
):
    result_path = pathlib.Path(results).resolve()
    fitname = result_path.name

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
        if at_least_nrep:
            log.warning(f"Postfit aiming for at least {nrep} replicas")
        else:
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
        passing_paths = filter_replicas(
            postfit_path, nnfit_path, fitname, chi2_threshold, arclength_threshold, integ_threshold
        )
        if len(passing_paths) < nrep:
            raise PostfitError("Number of requested replicas is too large")
        # Select the first nrep passing replicas
        if at_least_nrep:
            selected_paths = passing_paths
        else:
            selected_paths = passing_paths[:nrep]

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
        lhapdf.pathsPrepend(str(postfit_path))
        generatingPDF = PDF(fitname)
        lhio.generate_replica0(generatingPDF)

        # Test replica 0
        try:
            lhapdf.mkPDF(fitname, 0)
        except RuntimeError as e:
            raise PostfitError("CRITICAL ERROR: Failure in reading replica zero") from e
    log.info("\n\n*****************************************************************\n")
    log.info("Postfit complete")
    log.info("Please upload your results with:")
    log.info(f"\tvp-upload {result_path}\n")
    log.info("and install with:")
    log.info(f"\tvp-get fit {fitname}\n")
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
        '--at-least-nrep',
        action='store_true',
        help="nrep becomes the minimum number of required replicas. If there are more than nrep "
        "good replicas, all good replicas are written to the postfit folder.",
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
            args.at_least_nrep,
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
