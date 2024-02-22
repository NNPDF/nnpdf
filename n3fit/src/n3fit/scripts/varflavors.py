"""
varflavors.py

When producing a PDF with different maximal flavor number than the nf during the fit, this script
should be used. It copies the original fit folder after which it run ``evolven3fit``  on this fit
using the manually selected theory_id.

After running evolven3fit, this script replaces the ``AlphaS_MZ'' and ``MZ'' values in the .info
file, with the ``alphas`` and  ``Qref`` values from the theory database.

This script does not run postfit, that should still be done manually.
"""

from argparse import ArgumentParser
import logging
from pathlib import Path
import shutil
import subprocess
import sys

from validphys.loader import Loader
from validphys.renametools import change_name
from validphys.theorydbutils import fetch_theory

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

ll = Loader()
path_db = ll.theorydb_file


def main():
    parser = ArgumentParser(
        description="varflavors - a script to produce PDFs with flavor-number variations"
    )
    parser.add_argument(
        "fit_folder", help="Path to the folder containing the (pre-DGLAP) fit result"
    )
    parser.add_argument(
        "max_replicas",
        help="Maximum number of replicas on which to perform DGLAP evolution",
        type=int,
    )
    parser.add_argument("theory_id", help="ID of the theory used to evolve the fit", type=int)
    parser.add_argument(
        "--new_name",
        help="Create a copy of the input fit with this name before varying flavors of this new fit",
        type=str,
    )
    args = parser.parse_args()

    path_input_fit = Path(args.fit_folder)
    input_fit_name = path_input_fit.name
    path_info_file = path_input_fit / f"nnfit/{input_fit_name}.info"

    # 1. If the new_name flag is used, create a copy of the input fit
    if args.new_name:
        import tempfile

        path_output_fit = path_input_fit.with_name(args.new_name)
        if path_output_fit.exists():
            log.error("Destination path %s already exists.", path_output_fit.absolute())
            sys.exit(1)
        with tempfile.TemporaryDirectory(dir=path_input_fit.parent) as tmp:
            tmp = Path(tmp)
            copied_fit = tmp / input_fit_name
            shutil.copytree(path_input_fit, copied_fit, symlinks=True)
            newpath = change_name(copied_fit, args.new_name)
            newpath.rename(path_output_fit)
        log.info("Renaming with copy completed")
    else:
        if path_input_fit.exists():
            path_output_fit = path_input_fit
        else:
            log.error("Could not find fit. Path %s is not a directory.", path_input_fit.absolute())
            sys.exit(1)

    # 2. Run evolven3fit
    evolven3fit_command = (
        f"evolven3fit {path_output_fit} {args.max_replicas} --theory_id {args.theory_id} "
    )
    subprocess.run(evolven3fit_command, check=True, shell=True)
    log.info("DGLAP evolution completed")

    # 3. Overwrite the MZ and AlphaS_MZ values in the .info file with Qref and alphas values from
    # the theory db, repsectively
    theory_dict = fetch_theory(path_db, args.theory_id)
    Qref = theory_dict["Qref"]
    alphas = theory_dict["alphas"]

    path_temp_info = path_info_file.parent / "temp.info"
    with open(path_info_file, "r+") as f:
        with open(path_temp_info, "w") as tempfile:
            for line in f:
                if line.startswith("AlphaS_MZ:"):
                    line = f"AlphaS_MZ: {alphas}\n"
                if line.startswith("MZ"):
                    line = f"MZ: {Qref}\n"
                tempfile.write(line)
    path_temp_info.rename(path_info_file)
    log.info(
        "AlphaS_MZ and MZ in the .info file replace with alphas and Qref"
        "corresponding to theory %s in the theory db",
        args.theory_id,
    )
