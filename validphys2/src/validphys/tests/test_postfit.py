"""
test_postfit.py

Module for testing postfit.
"""
import json
from os import listdir

from validphys.scripts.postfit import main as postfit
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT
from reportengine.compat import yaml


def test_postfit():
    """Checks that the following happens when postfit is run on a pre-existing fit:
    - The postfit directory is generated
    - The expected files are generated in the PDF set
    - The replicas in the PDF set correspond correctly to the replicas in the fit
    - The number of PDF members is written to the info file correctly
    - The chi2 and arclength thresholds are correctly written to veto_count.json

    """
    # Load the fit and find its path
    l = Loader()
    fit = l.check_fit(FIT)
    fitpath = fit.path

    # Define arguments with which to run postfit afresh
    # Set thresholds to non-default values
    nrep = 2
    chi2_threshold = 3
    arclength_threshold = 5.2
    args = [
        str(nrep),
        str(fitpath),
        "--chi2-threshold",
        str(chi2_threshold),
        "--arclength-threshold",
        str(arclength_threshold),
    ]

    # Run postfit
    postfit(args)

    # Check that postfit directory is created
    postfitpath = fitpath / "postfit"
    assert postfitpath.is_dir()

    # Check that there are the expected files in the PDF set folder
    pdfsetpath = postfitpath / f"{FIT}"
    # Create set of expected files, inc. info file
    # Use sets so that the files are automatically sorted
    # Add one to nrep to account for replica 0
    expected_pdf_files = {f"{FIT}_{x:04d}.dat" for x in range(nrep + 1)}
    expected_pdf_files.add(f"{FIT}.info")
    # Find set of files that postfit actually generates
    generated_pdf_files = set(listdir(pdfsetpath))
    assert expected_pdf_files == generated_pdf_files

    # Check that replicas in the PDF set correspond to the fit replicas correctly
    for x in range(1, nrep):
        repnos = set()
        # [File in PDF set, file in fit]
        files = [pdfsetpath / f"{FIT}_{x:04d}.dat", postfitpath / f"replica_{x}/{FIT}.dat"]
        for file in files:
            with open(file, "r") as f:
                data = yaml.safe_load_all(f)
                metadata = next(data)
                repnos.add(metadata["FromMCReplica"])
        assert len(repnos) == 1

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
