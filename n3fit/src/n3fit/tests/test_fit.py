"""
    Regression tests for n3fit

    This file will run a fit with a runcard which includes:
        - A DIS dataset
        - A Hadronic dataset
        - Two positivity sets
    And checks that the results have not changed from the previous iteration of the code

    If the results are known to need a change,
    it is necessary to flag _something_ to regenerate regression
"""

import tempfile
import pathlib
import logging
import shutil
import os
import subprocess as sp

log = logging.getLogger(__name__)
REGRESSION_FOLDER = pathlib.Path().absolute() / "regressions"
QUICKNAME = "quickcard"
QUICKCARD = pathlib.Path().absolute() / f"{QUICKNAME}.yml"
EXE = "n3fit"
REPLICA = "1"


def load_data(path):
    """ Loads the info file of the fit into a list """
    info_file = path / f"{QUICKNAME}/nnfit/replica_{REPLICA}/{QUICKNAME}.fitinfo"
    f = open(info_file, "r")
    info = f.readlines()
    f.close()
    return info


def compare_lines(set1, set2):
    """ Returns true if set1 and set2 are equal """
    return set1 == set2


def test_fit():
    # create a /tmp folder
    tmp_name = tempfile.mkdtemp(prefix="nnpdf-")
    tmp_path = pathlib.Path(tmp_name)
    # cp runcard to tmp folder
    shutil.copy(QUICKCARD, tmp_path)
    os.chdir(tmp_path)
    # run the fit
    run_command = [EXE, QUICKCARD, REPLICA]
    sp.run(run_command)
    # read up the .dat files
    new_fitinfo = load_data(tmp_path)
    old_fitinfo = load_data(REGRESSION_FOLDER)
    # compare to the previous .dat file
    assert compare_lines(new_fitinfo, old_fitinfo)
