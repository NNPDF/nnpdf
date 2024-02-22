"""
vp-nextfitruncard

Command line tool to produce the runcard for an iterated fit given the input fit.

The script:
- Takes the input fit as a required argument and loads its runcard
- Updates the description of the fit interactively
- Uses the input fit as the t0 set
- Modifies the random seeds to values between 0 and 2**32 - 1 (max size of unsigned long int)
- Updates the preprocessing exponents
- Writes the runcard for the iterated fit to the current working directory, unless a different path
  is given as an argument
- Note that the runcard is written as long as it does not already exist in the path. This can be
  overridden by using the --force flag
"""

import argparse
import logging
import os
import pathlib
import sys

import prompt_toolkit

from reportengine import colors
from reportengine.compat import yaml
from validphys.api import API

# arguments for np.clip to enforce integrability.
# key should be identical to runcard key, first inner dictionary can contain
# either smallx or largex, the innermost dictionaries must be valid arguments
# for np.clip, this means BOTH a_min and a_max must be specified (even if one
# is left as None, indicating it is unconstrained.)
PREPROCESSING_LIMS = {
    "v": {"smallx": {"a_min": None, "a_max": 1.0}},
    "v3": {"smallx": {"a_min": None, "a_max": 1.0}},
    "v8": {"smallx": {"a_min": None, "a_max": 1.0}},
    "t3": {"smallx": {"a_min": None, "a_max": 1.0}},
    "t8": {"smallx": {"a_min": None, "a_max": 1.0}},
}


# Take command line arguments
def process_args():
    parser = argparse.ArgumentParser(description="Script to generate iterated fit runcard.")
    parser.add_argument("input_fit", help="Name of input fit.")
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=os.getcwd(),
        help="Directory to which the new runcard will be written. This must be a valid path. The default is the current working directory.",
    )
    parser.add_argument(
        "-f",
        "--force",
        help="If the runcard for the iterated fit already exists in the path, overwrite it.",
        action="store_true",
    )
    parser.add_argument(
        "--no-preproc-lims",
        action="store_true",
        help=(
            "Do not enforce any preprocessing constraints, which are chosen to "
            "ensure integrability. By default the following constraints are "
            f"used: {PREPROCESSING_LIMS}"
        ),
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
    # Convert given output directory to path and check it exists
    output_path = pathlib.Path(output_dir)
    if not output_path.is_dir():
        log.error("The specified output directory is not a valid path.")
        sys.exit(1)

    force = args.force
    if force:
        log.warning(
            "--force set to True. If the runcard for the iterated fit already exists in path to be "
            "written to, it will be overwritten."
        )

    output_fit = input_fit + "_iterated.yaml"
    runcard_path_out = output_path / output_fit
    # Check whether runcard with same name already exists in the path
    if runcard_path_out.exists() and not force:
        log.error(
            "Destination path %s already exists. If you wish to "
            "overwrite it, use the --force option.",
            runcard_path_out.absolute(),
        )
        sys.exit(1)

    # Update description of fit interactively
    description = API.fit(fit=input_fit).as_input()["description"]

    if not args.no_preproc_lims:
        preproc_lims = PREPROCESSING_LIMS
        log.info(
            "The following constraints will be used for preprocessing ranges, \n%s",
            yaml.dump(preproc_lims),
        )
    else:
        # don't enforce any limits.
        preproc_lims = None

    updated_description = interactive_description(description)

    iterated_runcard_yaml = API.iterated_runcard_yaml(
        fit=input_fit, _updated_description=updated_description, _flmap_np_clip_arg=preproc_lims
    )

    # Write new runcard to file
    with open(runcard_path_out, "w") as outfile:
        outfile.write(iterated_runcard_yaml)
        log.info("Runcard for iterated fit written to %s.", runcard_path_out.absolute())

    # Open new runcard with default editor, or if one is not set, with vi
    EDITOR = os.environ.get("EDITOR") if os.environ.get("EDITOR") else "vi"
    os.system(f"{EDITOR} {runcard_path_out}")


if __name__ == "__main__":
    main()
