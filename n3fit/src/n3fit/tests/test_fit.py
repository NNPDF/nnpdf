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

# this is needed for Travis to pass the test in mac
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
# make sure the GPU is not used
os.environ["CUDA_VISIBLE_DEVICES"] = ""
# https://keras.io/getting_started/faq/#how-can-i-obtain-reproducible-results-using-keras-during-development
os.environ["PYTHONHASHSEED"] = "0"

import pytest
import shutil
import pathlib
import logging
import subprocess as sp
from collections import namedtuple
from numpy.testing import assert_almost_equal
from reportengine.compat import yaml
import n3fit
from n3fit.performfit import initialize_seeds

log = logging.getLogger(__name__)
REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
EXE = "n3fit"
REPLICA = "1"
EXPECTED_MAX_FITTIME = 130  # seen mac ~ 180  and linux ~ 90


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
    """Returns true if the lines within set1 and set2 are the same
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


def auxiliary_performfit(tmp_path, timing=True):
    quickcard = f"{QUICKNAME}.yml"
    # Prepare the runcard
    quickpath = REGRESSION_FOLDER / quickcard
    weightpath = REGRESSION_FOLDER / "weights.h5"
    # read up the old info file
    old_fitinfo = load_data(REGRESSION_FOLDER / f"{QUICKNAME}.fitinfo")
    # cp runcard and weights to tmp folder
    shutil.copy(quickpath, tmp_path)
    shutil.copy(weightpath, tmp_path)
    # run the fit
    sp.run(f"{EXE} {quickcard} {REPLICA}".split(), cwd=tmp_path, check=True)
    # read up the .fitinfo files
    full_path = tmp_path / f"{QUICKNAME}/nnfit/replica_{REPLICA}/{QUICKNAME}.fitinfo"
    new_fitinfo = load_data(full_path)
    # compare to the previous .fitinfo file
    compare_lines(new_fitinfo[:5], old_fitinfo[:5], precision=1)
    # check that the times didnt grow in a weird manner
    time_path = tmp_path / f"{QUICKNAME}/nnfit/replica_{REPLICA}/{QUICKNAME}.time"
    if timing:
        # Better to catch up errors even when they happen to grow larger by chance
        with open(time_path, "r") as stream:
            times = yaml.safe_load(stream)
        fitting_time = times["walltime"]["replica_set_to_replica_fitted"]
        assert fitting_time < EXPECTED_MAX_FITTIME
    version_path = tmp_path / f"{QUICKNAME}/nnfit/replica_{REPLICA}/version.info"
    with open(version_path, "r") as stream:
        versions = yaml.safe_load(stream)
        assert versions["nnpdf"] == n3fit.__version__


@pytest.mark.darwin
def test_performfit(tmp_path):
    auxiliary_performfit(tmp_path, timing=False)


@pytest.mark.linux
def test_performfit_and_timing(tmp_path):
    auxiliary_performfit(tmp_path, timing=True)


def test_hyperopt(tmp_path):
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    # cp runcard to tmp folder
    shutil.copy(quickpath, tmp_path)
    # We just want to ensure that the hyperopt can run, but we need to kill it ourselves
    # 60 seconds should be enough
    with pytest.raises(sp.TimeoutExpired):
        sp.run(
            f"{EXE} {quickcard} {REPLICA} --hyperopt 1000".split(),
            cwd=tmp_path,
            timeout=60,
        )


def test_novalidation(tmp_path, timing=30):
    """ Runs a runcard without validation, success is assumed if it doesn't crash in 30 seconds """
    quickcard = f"noval-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    shutil.copy(quickpath, tmp_path)
    with pytest.raises(sp.TimeoutExpired):
        sp.run(f"{EXE} {quickcard} {REPLICA}".split(), cwd=tmp_path, timeout=timing)
