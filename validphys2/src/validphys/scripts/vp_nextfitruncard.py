"""
vp-nextfitruncard

Command line tool to produce the runcard for an iterated fit given the input fit.

The script:
- Takes the input fit as a required argument and loads its runcard
- Updated the description of the fit interactively
- Uses the input fit as the t0 set
- Modifies the random seeds to values between 0 and 1e10
- Writes the runcard for the iterated fit to the current working directory, unless a different path
  is given as an argument
- Note that the runcard is written as long as it does not already exist in the path. This can be
  overridden by using the --force flag
"""

import argparse
import os
import pathlib
import sys
import random
import logging
import prompt_toolkit

import NNPDF as nnpath

from reportengine.compat import yaml
from reportengine import colors
from reportengine.floatformatting import significant_digits

from validphys.eff_exponents import next_effective_exponents_table
from validphys.loader import Loader, PDFNotFound
from validphys.pdfbases import check_basis


# Take command line arguments
def process_args():
    parser = argparse.ArgumentParser(
        description="Script to generate iterated fit runcard."
    )
    parser.add_argument("input_fit", help="Name of input fit.")
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=os.getcwd(),
        help="Directory to which the new runcard will be written. The default is the current working directory.",
    )
    parser.add_argument(
        "-f",
        "--force",
        help="If the runcard for the iterated fit already exists in the path, overwrite it.",
        action="store_true",
    )
    args = parser.parse_args()
    return args


def interactive_description(original_description):
    """Set description of fit interactively. The description of the input fit is used as default."""
    default = original_description
    new_description = prompt_toolkit.prompt(
        "Enter a description for the new fit, taking into account that it should state that this fit is an iteration: \n",
        default=default,
    )
    if not new_description:
        return default
    return new_description


def main():
    # Logger for writing to screen
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.addHandler(colors.ColorHandler())

    args = process_args()

    input_fit = args.input_fit
    output_dir = args.output_dir
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
            "`vp-get fit <fit_name>`."
        )
        sys.exit(1)

    runcard_path_in = fit_path / "filter.yml"

    with open(runcard_path_in, "r") as infile:
        log.info(f"Input fit is {input_fit}.")
        runcard_data = yaml.load(infile, Loader=yaml.RoundTripLoader)
        log.info(f"Input runcard successfully read from {runcard_path_in.absolute()}.")

    # Update runcard with settings needed for iteration

    # Update description of fit interactively
    description = runcard_data["description"]
    updated_description = interactive_description(description)
    runcard_data["description"] = updated_description

    # Iterate t0
    runcard_data["datacuts"]["t0pdfset"] = input_fit

    # Update seeds with pseudorandom numbers between 0 and 1e10
    # Check if seeds exist especially since extra seeds needed in n3fit vs nnfit
    # Start with seeds in "fitting" section of runcard
    fitting_data = runcard_data["fitting"]
    fitting_seeds = ["seed", "trvlseed", "nnseed", "mcseed"]

    for seed in fitting_seeds:
        if seed in fitting_data:
            fitting_data[seed] = random.randrange(0, 1e10)

    # Next "closuretest" section of runcard
    closuretest_data = runcard_data["closuretest"]
    if "filterseed" in closuretest_data:
        closuretest_data["filterseed"] = random.randrange(0, 1e10)

    # Update preprocessing exponents
    # Firstly, find new exponents from PDF set that corresponds to fit
    basis = runcard_data["fitting"]["fitbasis"]
    checked = check_basis(basis, None)
    basis = checked["basis"]
    flavours = checked["flavours"]
    l = Loader()
    try:
        pdf = l.check_pdf(input_fit)
    except PDFNotFound:
        log.error(
            f"Could not find the PDF set {input_fit}. If the requested PDF set does exist, you can "
            "download it with `vp-get pdf <pdf_name>`."
        )
        sys.exit(1)
    new_exponents = next_effective_exponents_table(
        pdf=pdf, basis=basis, flavours=flavours
    )

    # Define function that we will use to round exponents to 4 significant figures
    rounder = lambda a: float(significant_digits(a, 4))

    # Update previous_exponents with new values
    previous_exponents = runcard_data["fitting"]["basis"]
    runcard_flavours = basis.to_known_elements(
        [ref_fl["fl"] for ref_fl in previous_exponents]
    ).tolist()
    for fl in flavours:
        alphas = new_exponents.loc[(f"${fl}$", r"$\alpha$")].values
        betas = new_exponents.loc[(f"${fl}$", r"$\beta$")].values
        previous_exponents[runcard_flavours.index(fl)]["smallx"] = [
            rounder(alpha) for alpha in alphas
        ]
        previous_exponents[runcard_flavours.index(fl)]["largex"] = [
            rounder(beta) for beta in betas
        ]

    # Start process of writing new runcard to file
    output_path = pathlib.Path(output_dir)
    output_fit = input_fit + "_iterated.yaml"
    runcard_path_out = output_path / output_fit

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

    # Open new runcard with default editor, or if one is not set, with vim
    EDITOR = os.environ.get("EDITOR") if os.environ.get("EDITOR") else "vim"
    os.system(f"{EDITOR} {runcard_path_out}")


if __name__ == "__main__":
    main()
