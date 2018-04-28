#!/usr/bin/env python
"""
    setup-fit - prepare and apply data cuts before fit (filter replacement)

    setup-fit constructs the fit [results] folder where data used by nnfit
    will be stored.
"""
import os
import sys
import re
import shutil
import argparse
import pathlib
import logging
import hashlib

from reportengine import colors

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(colors.ColorHandler())

RUNCARD_COPY_FILENAME = "filter.yml"
FILTER_OUTPUT_FOLDER = "filter"
MD5_FILENAME = "md5"


class SetupFitError(Exception):
    """Exception raised when setup-fit cannot succeed and knows why"""
    pass


def build_results_folder(runcard: str):
    """Parse runcard and create results folder"""
    file = pathlib.Path(runcard)
    # check file exists, is a file, has extension.
    if not file.exists():
        raise SetupFitError("Invalid runcard. File not found.")
    else:
        if not re.fullmatch(r'[\w.\-]+', runcard):
            raise SetupFitError("Invalid runcard. Must be alphanumeric.")
        if not file.is_file():
            raise SetupFitError("Invalid runcard. Must be a file.")

    # check filename
    filename, extension = os.path.splitext(runcard)
    if not len(extension):
        raise SetupFitError("Invalid runcard. File extension missing.")

    # check if results folder exists
    folder = pathlib.Path(filename)
    if folder.is_dir() and folder.exists():
        log.warning("Output folder already exists!")

    # create output folder
    output_folder = pathlib.Path(folder/FILTER_OUTPUT_FOLDER)
    output_folder.mkdir(mode=0o755, parents=True, exist_ok=True)

    # place a copy of runcard
    shutil.copy(runcard, folder/RUNCARD_COPY_FILENAME)
    return folder, output_folder


def store_md5(input_filename, output_filename):
    """Generate md5 key from file"""
    with open(input_filename, 'rb') as f:
        hash_md5 = hashlib.md5(f.read()).hexdigest()
    with open(output_filename, 'w') as g:
        g.write(hash_md5)
    log.info("md5 %s stored in %s" % (hash_md5, str(output_filename)))


def parse_command_line_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runcard', help="configuration file name")
    parser.add_argument('-d', '--debug', action='store_true', help='show debug messages')
    args = parser.parse_args()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    return args


def main():
    # parse command line arguments
    args = parse_command_line_args()

    try:
        # prepare results folder
        output_folder, filter_folder = build_results_folder(args.runcard)

        # load runcard

        # check for pdf sets

        # set rng seed for fake data

        # create md5 key
        store_md5(args.runcard, output_folder/MD5_FILENAME)

    except SetupFitError as e:
        log.error(f"Error in setup-fit:\n{e}")
        sys.exit(1)
    except Exception as e:
        log.critical(f"Bug in setup-fit ocurred. Please report it.")
        raise

    log.info("setup-fit completed successfully")
