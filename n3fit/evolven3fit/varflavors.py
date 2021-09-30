"""
varflavors.py

When producing a PDF with different maximal flavor number than the nf during the fit, we first run
evolven3fit by manually selecting a theory id of the theory correspodning to the desired nf
evolution.

After running evolven3fit, but before running posfit, this script should be run to replace the
``AlphaS_MZ'' and ``MZ'' values in the .info file, with the ``alphas'' and  ``Qref'' valuse from the
theory database.
"""

from argparse import ArgumentParser
from pathlib import Path

from validphys.loader import Loader
from validphys.theorydbutils import fetch_theory


ll = Loader()
path_db = ll.theorydb_file

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--fit_folder", required=True, help=f"Path to the folder containing the evolved fit"
    )
    parser.add_argument(
        "--theory_id", required=True, help=f"ID of the theory used to evolve the fit"
    )
    args = parser.parse_args()

    path_fit_folder = Path(args.fit_folder)
    fit_name = path_fit_folder.name
    path_info_file = path_fit_folder / f"nnfit/{fit_name}.info"

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
