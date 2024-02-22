"""
test_postfit.py

Module for testing postfit.
"""
import json
import os
import shutil
import subprocess as sp

from reportengine.compat import yaml
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT


def test_postfit(tmp):
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

    # Copy fit to temporary location and rename it to avoid any unwanted mutation of the PDF set in
    # the LHAPDF folder. Otherwise other tests may fail, for example
    shutil.copytree(fitpath, tmp / FIT)
    TMPFIT = "TEST"
    sp.run(f"vp-fitrename -c {fit.name} {TMPFIT}".split(), cwd=tmp, check=True)
    # Update fitpath so that it is the path of the fit that we will run postfit on
    fitpath = tmp / TMPFIT

    # Define arguments with which to run postfit afresh
    # Set thresholds to non-default values
    nrep = 2
    chi2_threshold = 3
    arclength_threshold = 5.2
    integrability_threshold = 1.0

    # Run postfit
    sp.run(
        f"postfit {nrep} {TMPFIT} --chi2-threshold {chi2_threshold} --arclength-threshold {arclength_threshold} --integrability-threshold {integrability_threshold}".split(),
        cwd=tmp,
        check=True,
    )

    # Check that postfit directory is created
    postfitpath = fitpath / "postfit"
    assert (
        postfitpath.is_dir()
    ), f"The postfit directory has not been created as expected at {postfitpath}."

    # Check that there are the expected files in the PDF set folder
    pdfsetpath = postfitpath / f"{TMPFIT}"
    # Create set of expected files, inc. info file
    # Use sets so that the files are automatically sorted
    # Start counting from zero because of replica 0
    # Add one to nrep to account for the last replica
    expected_pdf_files = {f"{TMPFIT}_{x:04d}.dat" for x in range(nrep + 1)}
    expected_pdf_files.add(f"{TMPFIT}.info")
    # Find set of files that postfit actually generates
    generated_pdf_files = set(os.listdir(pdfsetpath))
    assert (
        expected_pdf_files == generated_pdf_files
    ), f"""The set of files generated for the PDF set by postfit differs from the set of expected files.
           The problematic files are {", ".join(expected_pdf_files.symmetric_difference(generated_pdf_files))})."""

    # Check that replicas in the PDF set correspond to the fit replicas correctly
    for x in range(1, nrep + 1):
        repnos = set()
        # [File in PDF set, file in fit]
        files = [pdfsetpath / f"{TMPFIT}_{x:04d}.dat", postfitpath / f"replica_{x}/{TMPFIT}.dat"]
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
    infopath = postfitpath / f"{TMPFIT}/{TMPFIT}.info"
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
        assert (
            veto_count["integrability_threshold"] == integrability_threshold
        ), f"Postfit has not written the integrability threshold correctly to {vetopath}."
