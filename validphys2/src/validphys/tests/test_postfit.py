"""
test_postfit.py

Module for testing postfit.
"""
import json
from os import listdir

from validphys.scripts.postfit import main as postfit
from validphys.loader import FallbackLoader as Loader
from reportengine.compat import yaml

# Use fit that is only used by this test, rather than importing the default from conftest.py
# We do this to avoid any unwanted interference between the tests, in particular because this
# test modifies the postfit folder of the fit and therefore the info associated with its LHAPDF set
FIT = "191015-mw-001_for_postfit_test"

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
    assert (
        postfitpath.is_dir()
    ), f"The postfit directory has not been created as expected at {postfitpath}."

    # Check that there are the expected files in the PDF set folder
    pdfsetpath = postfitpath / f"{FIT}"
    # Create set of expected files, inc. info file
    # Use sets so that the files are automatically sorted
    # Start counting from zero because of replica 0
    # Add one to nrep to account for the last replica
    expected_pdf_files = {f"{FIT}_{x:04d}.dat" for x in range(nrep + 1)}
    expected_pdf_files.add(f"{FIT}.info")
    # Find set of files that postfit actually generates
    generated_pdf_files = set(listdir(pdfsetpath))
    assert (
        expected_pdf_files == generated_pdf_files
    ), f"""The set of files generated for the PDF set by postfit differs from the set of expected files.
           The problematic files are {", ".join(expected_pdf_files.symmetric_difference(generated_pdf_files))})."""

    # Check that replicas in the PDF set correspond to the fit replicas correctly
    for x in range(1, nrep + 1):
        repnos = set()
        # [File in PDF set, file in fit]
        files = [pdfsetpath / f"{FIT}_{x:04d}.dat", postfitpath / f"replica_{x}/{FIT}.dat"]
        for file in files:
            with open(file, "r") as f:
                data = yaml.safe_load_all(f)
                metadata = next(data)
                repnos.add(metadata["FromMCReplica"])
        assert (
            len(repnos) == 1
        ), f"""There is a mismatch between the replica number written to the PDF set by postfit and
               the fit replica it recorded for replica_{x}."""

    # Check that number of PDF members is written correctly
    infopath = postfitpath / f"{FIT}/{FIT}.info"
    with open(infopath, "r") as f:
        data = yaml.safe_load(f)
        # Add one to nrep to account for replica 0
        assert (
            data["NumMembers"] == nrep + 1
        ), f"Postfit has not written the number of PDF members correctly to {infopath}."

    # Check that chi2 and arclength thresholds are recorded correctly
    vetopath = postfitpath / "veto_count.json"
    with open(vetopath, "r") as f:
        veto_count = json.load(f)
        assert (
            veto_count["chi2_threshold"] == chi2_threshold
        ), f"Postfit has not written the chi2 threshold correctly to {vetopath}."
        assert (
            veto_count["arclength_threshold"] == arclength_threshold
        ), f"Postfit has not written the arclength threshold correctly to {vetopath}."
