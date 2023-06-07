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

import json
import pytest
import shutil
import pathlib
import logging
import subprocess as sp
from numpy.testing import assert_equal, assert_allclose
import n3fit
from validphys.n3fit_data import replica_trvlseed, replica_nnseed, replica_mcseed

log = logging.getLogger(__name__)
REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
EXE = "n3fit"
REPLICA = "1"
EXPECTED_MAX_FITTIME = 130  # seen mac ~ 180  and linux ~ 90


def load_data(info_file):
    """Loads the information of the fit from the json files"""
    return json.loads(info_file.read_text(encoding="utf-8"))


def test_initialize_seeds():
    # Regression tests for seed generation
    replicas = [1, 4]

    trvlseeds = [replica_trvlseed(rep, 4) for rep in replicas]
    assert trvlseeds == [2005877882, 741720773]
    nnseeds = [replica_nnseed(rep, 7) for rep in replicas]
    assert nnseeds == [327741615, 1369975286]
    mcseeds = [replica_mcseed(rep, 1, True) for rep in replicas]
    assert mcseeds == [1791095845, 1857819720]

    mcseeds_nomc = [replica_mcseed(rep, 1000, False) for rep in replicas]
    assert mcseeds_nomc == [None, None]

    # test that we always get same answer
    same_replicas = [3]*10
    assert len({replica_trvlseed(rep, 4) for rep in same_replicas}) == 1
    assert len({replica_nnseed(rep, 7) for rep in same_replicas}) == 1
    assert len({replica_mcseed(rep, 1, True) for rep in same_replicas}) == 1


def auxiliary_performfit(tmp_path, replica=1, timing=True, rel_error=2e-3):
    """Fits quickcard and checks the json file to ensure the results have not changed.
    """
    quickcard = f"{QUICKNAME}.yml"
    # Prepare the runcard
    quickpath = REGRESSION_FOLDER / quickcard
    weightpath = REGRESSION_FOLDER / f"weights_{replica}.h5"
    # read up the previous json file for the given replica
    old_json = load_data(REGRESSION_FOLDER / f"{QUICKNAME}_{replica}.json")
    # cp runcard and weights to tmp folder
    shutil.copy(quickpath, tmp_path)
    shutil.copy(weightpath, tmp_path / "weights.h5")
    # run the fit
    sp.run(f"{EXE} {quickcard} {replica}".split(), cwd=tmp_path, check=True)
    # read up json files
    full_json = tmp_path / f"{QUICKNAME}/nnfit/replica_{replica}/{QUICKNAME}.json"
    new_json = load_data(full_json)
    # Now compare to regression results, taking into account precision won't be 100%
    equal_checks = ["stop_epoch", "pos_state"]
    approx_checks = ["erf_tr", "erf_vl", "chi2", "best_epoch", "arc_lengths", "integrability", "best_epoch"]
    for key in equal_checks:
        assert_equal(new_json[key], old_json[key])
    for key in approx_checks:
        if old_json[key] is None and new_json[key] is None:
            continue
        assert_allclose(new_json[key], old_json[key], rtol=rel_error)
    # check that the times didnt grow in a weird manner
    if timing:
        # Better to catch up errors even when they happen to grow larger by chance
        times = new_json["timing"]
        fitting_time = times["walltime"]["replica_set_to_replica_fitted"]
        assert fitting_time < EXPECTED_MAX_FITTIME
    # For safety, check also the version
    assert new_json["version"]["nnpdf"] == n3fit.__version__


@pytest.mark.darwin
def test_performfit(tmp_path):
    auxiliary_performfit(tmp_path, replica=2, timing=False, rel_error=1e-1)


@pytest.mark.linux
@pytest.mark.parametrize("replica", [1, 2])
def test_performfit_and_timing(tmp_path, replica):
    auxiliary_performfit(tmp_path, replica=replica, timing=True)


@pytest.mark.skip(reason="Still not implemented in parallel mode")
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
    """ Runs a runcard without validation, success is assumed if it doesn't crash in 30 seconds 
    Checks that the code is able to run when the trvl frac is set to 1.0
    """
    quickcard = f"noval-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    shutil.copy(quickpath, tmp_path)
    with pytest.raises(sp.TimeoutExpired):
        sp.run(f"{EXE} {quickcard} {REPLICA}".split(), cwd=tmp_path, timeout=timing)


def test_weirdbasis(tmp_path, timing=30):
    """ Runs a runcard with perturbative charm basis but an intrinsic charm theory
    so the run will be stopped by the checks """
    # Once we have a light theory with perturbative charm for testing, this test can be enabled
    # to do the commented docstring
#     """ Runs a runcard with perturbative charm, success is assumed if it doesn't crash in 30 seconds
#     Checks that the code runs when it needs to rotate the output of the NN to the NN31IC basis
#     """
    quickcard = f"pc-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    shutil.copy(quickpath, tmp_path)
#     with pytest.raises(sp.TimeoutExpired):
    with pytest.raises(sp.CalledProcessError):
        sp.run(f"{EXE} {quickcard} {REPLICA}".split(), cwd=tmp_path, timeout=timing, check=True)
