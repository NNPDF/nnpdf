"""
test_postfit.py

Module for testing postfit.
"""
import json

from validphys.scripts.postfit import main as postfit
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT
from reportengine.compat import yaml


def test_postfit():
    """Checks that the following happens when postfit is run on a pre-existing fit:
    - The names of the files that postfit generates in a fresh run match those in the original run
    - The number of PDF members is written to the info file correctly
    - The chi2 and arclength thresholds are correctly written to veto_count.json

    """
    # Load the fit and find its path
    l = Loader()
    fit = l.check_fit(FIT)
    fitpath = fit.path

    # Record names of files in original run of postfit
    postfitpath = fitpath / "postfit"
    original_postfit_files = [x for x in postfitpath.glob("**/*") if x.is_file()]

    # Define arguments with which to run postfit afresh
    nrep = 30
    chi2_threshold = 3
    arclength_threshold = 5
    args = [
        str(nrep),
        str(fitpath),
        "--chi2-threshold",
        str(chi2_threshold),
        "--arclength-threshold",
        str(arclength_threshold),
    ]

    # Rerun postfit
    postfit(args)

    # Check that postfit directory is created
    assert postfitpath.is_dir()

    # Check that the files in the fresh run of postfit match those in the original run
    regenerated_postfit_files = [x for x in postfitpath.glob("**/*") if x.is_file()]
    assert original_postfit_files == regenerated_postfit_files

    # Check that number of PDF members is written correctly
    with open(postfitpath / f"{FIT}/{FIT}.info", "r") as f:
        data = yaml.safe_load(f)
        # Add one to nrep to account for replica 0
        assert data["NumMembers"] == nrep + 1

    # Check that chi2 and arclength thresholds are recorded correctly
    with open(postfitpath / f"veto_count.json", "r") as f:
        veto_count = json.load(f)
        assert veto_count["chi2_threshold"] == chi2_threshold
        assert veto_count["arclength_threshold"] == arclength_threshold
