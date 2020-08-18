"""
vp-nextfitruncard

Command line tool to produce the runcard for an iterated fit given the non-iterated fit.

The script:
- Takes the non-iterated fit as a required argument and loads its runcard
- Uses the non-iterated fit as the t0 set
- Modifies the random seed to a value between 0 and 1e10
- Writes the runcard for the iterated fit to the user's NNPDF config path as long as it does not
  already exist, unless the --force flag is given
"""

import argparse
import pathlib
import sys
import random
import logging

import NNPDF as nnpath

from reportengine.compat import yaml
from reportengine import colors


# Take command line arguments
def process_args():
    parser = argparse.ArgumentParser(
        description="Script to generate iterated fit runcard."
    )
    parser.add_argument("input_fit", help="Name of input fit.")
    parser.add_argument(
        "-f",
        "--force",
        help="If the runcard for the iterated fit already exists in the path, overwrite it.",
        action="store_true",
    )
    args = parser.parse_args()
    return args


def main():
    # Logger for writing to screen
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.addHandler(colors.ColorHandler())

    args = process_args()

    input_fit = args.input_fit
    force = args.force
    if force:
        log.warning(
            "--force set to True. If the runcard for the iterated fit already exists in config "
            "path, it will be overwritten."
        )

    results_path = pathlib.Path(nnpath.get_results_path())
    fit_path = results_path / input_fit

    if not fit_path.is_dir():
        log.error(
            "Could not find the specified fit. The following path is not a directory: "
            f"{fit_path.absolute()}. If the requested fit does exist, you can download it with "
            "`vp-get fit <fit_name>`"
        )
        sys.exit(1)

    runcard_path_in = fit_path / "filter.yml"

    with open(runcard_path_in, "r") as infile:
        log.info(f"Input fit is {input_fit}.")
        runcard_data = yaml.load(infile, Loader=yaml.RoundTripLoader)
        log.info(f"Input runcard successfully read from {runcard_path_in.absolute()}.")

    # Update runcard with settings needed for iteration
    # Iterate t0
    runcard_data["datacuts"]["t0pdfset"] = input_fit
    # Update seed with pseudorandom number between 0 and 1e10
    runcard_data["fitting"]["seed"] = random.randrange(0, 1e10)

    config_path = pathlib.Path(nnpath.get_config_path())
    output_fit = input_fit + "_iterated.yaml"
    runcard_path_out = config_path / output_fit

    if runcard_path_out.exists() and not force:
        log.error(
            f"Destination path {runcard_path_out.absolute()} already exists. If you wish to "
            "overwrite it, use the --force option."
        )
        sys.exit(1)

    with open(runcard_path_out, "w") as outfile:
        log.info("Dumping runcard for iterated fit.")
        yaml.dump(runcard_data, outfile, Dumper=yaml.RoundTripDumper)
        log.info(f"Runcard for iterated fit written to {runcard_path_out.absolute()}.")
