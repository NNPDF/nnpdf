#!/usr/bin/env python
"""
vp-pdffromreplicas

Take a pre-existing source ``PDF`` with MC replicas and create a new ``PDF``
with replicas subsampled from the source ``PDF``. Replicas will be sampled
uniformly.

source ``PDF`` will be downloaded if it cannot be found locally.

This script will not overwrite any existing files, so a ``PDF`` cannot already
exist with the same name as the output ``PDF``.

Note that whilst in principle it is possible to create a single replica ``PDF``
whose replica 0 is simply the same as replica 1, LHAPDF has a restriction which
requires ``PDF`` s to have at least 2 replicas (plus replica 0). To handle this
special case if ``replicas == 1``, then replica 2 will be a duplicate of replica
1, satisfying the minimum number of replicas whilst retaining the property
that replica 1 and replica 0 are identical.

"""

import argparse
import logging
import pathlib
import random
import shutil
import sys
import tempfile

import pandas as pd

from reportengine import colors
from reportengine.compat import yaml
from validphys import lhaindex
from validphys.lhio import new_pdf_from_indexes
from validphys.loader import FallbackLoader

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())


def check_none_or_gt_one(value):
    """Check the ``value`` supplied can be interpreted as an integer and is greater
    than one.

    Returns
    -------
    int
        supplied value cast to integer.
    """

    try:
        ivalue = int(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"{value} cannot be interpreted as an integer.") from e
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is an invalid positive int value.")
    return ivalue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_pdf",
        type=str,
        help="The name of the source PDF, from which the new PDF's replicas will be sampled.",
    )
    parser.add_argument(
        "replicas",
        type=check_none_or_gt_one,
        help="Number of replicas to sample, this number should not be greater than the number of replicas <input_pdf> has.",
    )
    parser.add_argument(
        "--output-name",
        "-o",
        type=str,
        default=None,
        help="Output name for the resulting PDF, defaults to <input_pdf>_<replicas>",
    )
    parser.add_argument(
        "--save-indices",
        "-s",
        action="store_true",
        help="Flag to save a CSV with a mapping of new replica indices to old replica indices.",
    )
    args = parser.parse_args()

    loader = FallbackLoader()
    input_pdf = loader.check_pdf(args.input_pdf)

    if input_pdf.error_type != "replicas":
        log.error("Error type of input PDF must be `replicas` not `%s`", input_pdf.error_type)
        sys.exit(1)

    if args.replicas > len(input_pdf) - 1:
        log.error(
            "Too many replicas requested: %s. The source PDF has %s",
            args.replicas,
            len(input_pdf) - 1,
        )
        sys.exit(1)

    indices = random.sample(range(1, len(input_pdf)), k=args.replicas)

    if args.output_name is None:
        output_name = args.input_pdf + f"_{args.replicas}"
    else:
        output_name = args.output_name

    with tempfile.TemporaryDirectory() as f:
        try:
            new_pdf_from_indexes(
                input_pdf, indices, set_name=output_name, folder=pathlib.Path(f), installgrid=True
            )
        except FileExistsError:
            log.error(
                "A PDF is already installed at %s, consider choosing a different output name.",
                pathlib.Path(lhaindex.get_lha_datapath()) / output_name,
            )
            sys.exit(1)

    if args.replicas == 1:
        log.warning(
            "PDFs in the LHAPDF format are required to have 2 replicas, copying "
            "replica 1 to replica 2"
        )
        base_name = str(pathlib.Path(lhaindex.get_lha_datapath()) / output_name / output_name)

        shutil.copyfile(base_name + "_0001.dat", base_name + "_0002.dat")
        # fixup info file
        with open(base_name + ".info", "r") as f:
            info_file = yaml.safe_load(f)

        info_file["NumMembers"] = 3
        with open(base_name + ".info", "w") as f:
            yaml.dump(info_file, f)

        # here we update old indices in case the user creates
        # the original_index_mapping.csv
        indices = 2 * indices

    if args.save_indices:
        index_file = (
            pathlib.Path(lhaindex.get_lha_datapath()) / output_name / "original_index_mapping.csv"
        )
        log.info("Saving output PDF/input PDF replica index mapping to %s", index_file)
        with open(index_file, "w+") as f:
            pd.DataFrame(
                list(enumerate(indices, 1)),
                columns=[f"{output_name} replica index", f"{args.input_pdf} replica index"],
            ).to_csv(f, index=False)


if __name__ == "__main__":
    main()
