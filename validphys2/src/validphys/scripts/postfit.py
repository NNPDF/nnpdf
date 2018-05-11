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

import sys
import os.path
import shutil
import pathlib
import argparse
import itertools
from glob import glob
import logging
import tempfile

import lhapdf

from reportengine import colors
from validphys import lhio
from validphys import fitdata
from validphys import fitveto
from validphys.core import PDF

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
    """ Forms a relative symlink between 'source' and 'dest' """
    relativepath = os.path.relpath(source, dest.parents[0])
    os.symlink(relativepath, dest)


def set_lhapdf_info(info_path, nrep):
    """ Sets the LHAPDF info file NumMembers field"""
    with open(info_path, 'r+') as f:
        txt = f.read()
        f.seek(0)
        f.write(txt.replace('REPLACE_NREP', str(nrep+1)))
        f.truncate()


class PostfitError(Exception):
    """Exception raised when postfit cannot suceed and knows why"""
    pass

class FatalPostfitError(Exception):
    """Excption raised when some corrupted input is detected"""
    pass

def filter_replicas(postfit_path, nnfit_path, fitname):
    """ Find the paths of all replicas passing the standard NNPDF fit vetoes
    as defined in fitveto.py. Returns a list of the replica directories that pass."""
    # This glob defines what is considered a valid replica
    # all the following code uses paths from this glob
    all_replicas   = glob(f"{nnfit_path}/replica_*/")
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
                f"Error when loading replica information:\n {e}") from e
    fit_vetoes = fitveto.determine_vetoes(fitinfo)
    fitveto.save_vetoes(fit_vetoes, postfit_path / 'veto_count.json')

    for key in fit_vetoes:
        log.info("%d replicas pass %s" % (sum(fit_vetoes[key]), key))
    passing_paths = list(itertools.compress(valid_paths, fit_vetoes["Total"]))
    return passing_paths


def postfit(results: str, nrep: int):
    result_path = pathlib.Path(results).resolve()
    fitname = result_path.name

    # Paths
    nnfit_path   = result_path / 'nnfit'    # Path of nnfit replica output
    # Create a temporary path to store work in progress and move it to
    # the final location in the end,
    postfit_path = pathlib.Path(tempfile.mkdtemp(prefix='postfit_work_deleteme_',
        dir=result_path))
    final_postfit_path = result_path / 'postfit'
    LHAPDF_path  = postfit_path/fitname     # Path for LHAPDF grid output

    if not fitdata.check_nnfit_results_path(result_path):
        raise PostfitError('Postfit cannot find a valid results path')
    if not fitdata.check_lhapdf_info(result_path, fitname):
        raise PostfitError('Postfit cannot find a valid LHAPDF info file')

    nrep = int(nrep)
    log.warn("Postfit aiming for %d replicas" % nrep)

    # Generate postfit and LHAPDF directory
    if final_postfit_path.is_dir():
        log.warn(f"WARNING: Removing existing postfit directory: {postfit_path}")
        shutil.rmtree(final_postfit_path)
    os.mkdir(LHAPDF_path)

    # Setup postfit log
    postfitlog = logging.FileHandler(postfit_path/'postfit.log', mode='w')
    log.addHandler(postfitlog)

    # Perform postfit selection
    passing_paths = filter_replicas(postfit_path, nnfit_path, fitname)
    if len(passing_paths) < nrep:
        raise PostfitError("Number of requested replicas is too large")
    # Select the first nrep passing replicas
    selected_paths = passing_paths[:nrep]


    # Copy info file
    info_source_path = nnfit_path.joinpath(f'{fitname}.info')
    info_target_path = LHAPDF_path.joinpath(f'{fitname}.info')
    shutil.copy2(info_source_path, info_target_path)
    set_lhapdf_info(info_target_path, nrep)

    # Generate symlinks
    for drep, source_path in enumerate(selected_paths, 1):
        # Symlink results to postfit directory
        source_dir = pathlib.Path(source_path).resolve()
        target_dir = postfit_path.joinpath('replica_%d' % drep)
        relative_symlink(source_dir, target_dir)
        # Symlink results to pdfset directory
        source_grid = source_dir.joinpath(fitname+'.dat')
        target_file = f'{fitname}_{drep:04d}.dat'
        target_grid = LHAPDF_path.joinpath(target_file)
        relative_symlink(source_grid, target_grid)

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
    postfit_path.rename(final_postfit_path)
    log.info("\n\n*****************************************************************\n")
    log.info("Postfit complete")
    log.info("Please upload your results with:")
    log.info(f"\tvp-upload --fit {result_path}\n")
    log.info("and install with:")
    log.info(f"\tvp-get fit {fitname}\n")
    log.info("*****************************************************************\n\n")


def main():
    splash()
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('nrep', type=int, help="Number of desired replicas")
    parser.add_argument('result_path', help="Folder containig the "
                                            "results of the fit")
    parser.add_argument('-d', '--debug', action='store_true', help='show debug messages')
    args = parser.parse_args()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    try:
        postfit(args.result_path, args.nrep)
    except PostfitError as e:
        log.error(f"Error in postfit:\n{e}")
        sys.exit(1)
    except FatalPostfitError as e:
        log.error(f"Corrupted input encountered")
        print(
            colors.color_exception(e.__class__, e, e.__traceback__),
            file=sys.stderr)
    except Exception as e:
        log.critical(f"Bug in postfit ocurred. Please report it.")
        print(
            colors.color_exception(e.__class__, e, e.__traceback__),
            file=sys.stderr)
