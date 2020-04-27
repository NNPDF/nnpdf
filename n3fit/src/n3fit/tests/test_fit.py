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

import os
import shutil
import pathlib
import logging
import tempfile
import subprocess as sp
from collections import namedtuple
from numpy.testing import assert_almost_equal
from n3fit.performfit import initialize_seeds

log = logging.getLogger(__name__)
REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
QUICKCARD = f"{QUICKNAME}.yml"
QUICKPATH = REGRESSION_FOLDER / QUICKCARD
EXE = "n3fit"
REPLICA = "1"


def load_data(info_file):
    """ Loads the info file of the fit into a list """
    with open(info_file, "r") as f:
        info = f.read()
        return info.split()


def compare_two(val1, val2, precision=6):
    """ Compares value 1 and value 2 attending to their type """
    try:
        num_1 = float(val1)
        num_2 = float(val2)
        assert_almost_equal(num_1, num_2, decimal=precision)
    except ValueError:
        assert val1 == val2


def compare_lines(set1, set2, precision=6):
    """ Returns true if the lines within set1 and set2 are the same
    The numbers are compared up to `precision`
    """
    for val1, val2 in zip(set1, set2):
        compare_two(val1, val2, precision=precision)


def test_initialize_seeds():
    # Regression tests for seed generation
    replicas = [1, 4]
    Seeds = namedtuple("Seeds", ["trvlseeds", "nnseeds", "mcseeds"])
    regression = Seeds(
        trvlseeds=[2005877882, 741720773],
        nnseeds=[327741615, 87982037],
        mcseeds=[1791095845, 2029572362],
    )
    result = initialize_seeds(replicas, 4, 7, 1, True)
    assert result == regression
    result_nomc = initialize_seeds(replicas, 10, 100, 1000, False)
    regression_nomc = Seeds(
        trvlseeds=[1165313289, 2124247567], nnseeds=[186422792, 1315999533], mcseeds=[]
    )
    assert result_nomc == regression_nomc


def test_performfit():
    return
    # read up the old info file
    old_fitinfo = load_data(REGRESSION_FOLDER / f"{QUICKNAME}.fitinfo")
    # create a /tmp folder
    tmp_name = tempfile.mkdtemp(prefix="nnpdf-")
    tmp_path = pathlib.Path(tmp_name)
    # cp runcard to tmp folder
    shutil.copy(QUICKPATH, tmp_path)
    # run the fit
    sp.run(f"{EXE} {QUICKCARD} {REPLICA}".split(), cwd=tmp_path, check=True)
    # read up the .fitinfo files
    full_path = tmp_path / f"{QUICKNAME}/nnfit/replica_{REPLICA}/{QUICKNAME}.fitinfo"
    new_fitinfo = load_data(full_path)
    # compare to the previous .fitinfo file
    compare_lines(new_fitinfo[:5], old_fitinfo[:5], precision=1)
