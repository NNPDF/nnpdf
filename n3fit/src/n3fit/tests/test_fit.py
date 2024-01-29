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

import json
import logging
import pathlib
import shutil
import subprocess as sp

import h5py
import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

import n3fit
from reportengine.compat import yaml
from validphys.n3fit_data import replica_mcseed, replica_nnseed, replica_trvlseed

log = logging.getLogger(__name__)
REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
QUICKNAME_QED = "quickcard_qed"
QUICKNAME_SEQUENTIAL = "quickcard-sequential"
QUICKNAME_PARALLEL = "quickcard-parallel"
EXE = "n3fit"
REPLICA = "1"
EXPECTED_MAX_FITTIME = 130  # seen mac ~ 180  and linux ~ 90


def _load_json(info_file):
    """Loads the information of the fit from the json files"""
    return json.loads(info_file.read_text(encoding="utf-8"))


def _load_exportgrid(exportgrid_file):
    """Loads the exportgrid file"""
    return yaml.safe_load(exportgrid_file.read_text())


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
    same_replicas = [3] * 10
    assert len({replica_trvlseed(rep, 4) for rep in same_replicas}) == 1
    assert len({replica_nnseed(rep, 7) for rep in same_replicas}) == 1
    assert len({replica_mcseed(rep, 1, True) for rep in same_replicas}) == 1


def check_fit_results(
    base_path, fitname, replica, regression_json, regenerate=False, rel_error=2e-3, timing=False
):
    """Regression test checker, checks that the given fit produces the right
    json and exportgrid files.

    This function will check that the fit ``fitname`` run in the directory ``base_path``
    produces the same values in the ``.json`` and ``.exportgrid`` files for replica ``replica``
    as can be read from the file ``regression_json``.
    The regression exportgrid is understood to be the same as the ``.json`` file with the
    extension changed to ``.exportgrid``.

    If ``regenerate`` is set to True, it will generate new files instead of testing
    """
    new_json_file = base_path / f"{fitname}/nnfit/replica_{replica}/{fitname}.json"

    new_expgrid_file = new_json_file.with_suffix(".exportgrid")
    old_expgrid_file = regression_json.with_suffix(".exportgrid")

    if regenerate:
        shutil.copy2(new_expgrid_file, old_expgrid_file)
        shutil.copy2(new_json_file, regression_json)

    # Compare json results
    new_json = _load_json(new_json_file)
    old_json = _load_json(regression_json)

    equal_checks = ["stop_epoch", "pos_state"]
    approx_checks = ["erf_tr", "erf_vl", "chi2", "best_epoch", "best_epoch"]
    relaxed_checks = ["arc_lengths", "integrability"]
    for key, value in new_json.items():
        reference = old_json[key]
        err_msg = f"error for .json: {key}"
        if key in equal_checks:
            assert_equal(value, reference, err_msg=err_msg)
        elif key in approx_checks:
            assert_allclose(value, reference, err_msg=err_msg, rtol=rel_error)
        elif key in relaxed_checks:
            assert_allclose(value, reference, err_msg=err_msg, rtol=rel_error * 10)
        elif key == "preprocessing":
            for ref, cur in zip(reference, value):
                err_msg += f" - {ref['fl']}"
                assert_allclose(ref["smallx"], cur["smallx"], err_msg=err_msg, rtol=rel_error)
                assert_allclose(ref["largex"], cur["largex"], err_msg=err_msg, rtol=rel_error)

    # check that the times didnt grow in a weird manner
    if timing:
        # Better to catch up errors even when they happen to grow larger by chance
        times = new_json["timing"]
        fitting_time = times["walltime"]["replica_set_to_replica_fitted"]
        assert fitting_time < EXPECTED_MAX_FITTIME

    # For safety, check also the version
    assert new_json["version"]["nnpdf"] == n3fit.__version__

    new_expgrid = _load_exportgrid(new_expgrid_file)
    old_expgrid = _load_exportgrid(old_expgrid_file)

    # Now compare the exportgrids
    for key, value in new_expgrid.items():
        reference = old_expgrid[key]
        err_msg = f"error for .exportgrid: {key}"
        if key == "pdfgrid":
            assert_allclose(value, reference, rtol=rel_error, atol=1e-6, err_msg=err_msg)
        else:
            assert_equal(value, reference, err_msg=err_msg)


def _auxiliary_performfit(tmp_path, runcard=QUICKNAME, replica=1, timing=True, rel_error=2e-3):
    """Fits quickcard and checks the json file to ensure the results have not changed."""
    quickcard = f"{runcard}.yml"
    # Prepare the runcard
    quickpath = REGRESSION_FOLDER / quickcard
    weightpath = REGRESSION_FOLDER / f"weights_{replica}.h5"
    # read up the previous json file for the given replica
    old_json_file = REGRESSION_FOLDER / f"{runcard}_{replica}.json"
    # cp runcard and weights to tmp folder
    shutil.copy(quickpath, tmp_path)
    shutil.copy(weightpath, tmp_path / "weights.h5")
    # run the fit
    sp.run(f"{EXE} {quickcard} {replica}".split(), cwd=tmp_path, check=True)

    # And compare
    check_fit_results(tmp_path, runcard, replica, old_json_file, timing=timing, rel_error=rel_error)


@pytest.mark.darwin
@pytest.mark.parametrize("runcard", [QUICKNAME, QUICKNAME_QED])
def test_performfit(tmp_path, runcard):
    _auxiliary_performfit(tmp_path, runcard=runcard, replica=3, timing=False, rel_error=1e-1)


@pytest.mark.linux
@pytest.mark.parametrize("replica", [1, 3])
@pytest.mark.parametrize("runcard", [QUICKNAME, QUICKNAME_QED])
def test_performfit_and_timing(tmp_path, runcard, replica):
    _auxiliary_performfit(tmp_path, runcard=runcard, replica=replica, timing=True)


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
        sp.run(f"{EXE} {quickcard} {REPLICA} --hyperopt 1000".split(), cwd=tmp_path, timeout=60)


def test_novalidation(tmp_path, timing=30):
    """Runs a runcard without validation, success is assumed if it doesn't crash in 30 seconds
    Checks that the code is able to run when the trvl frac is set to 1.0
    """
    quickcard = f"noval-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    shutil.copy(quickpath, tmp_path)
    with pytest.raises(sp.TimeoutExpired):
        sp.run(f"{EXE} {quickcard} {REPLICA}".split(), cwd=tmp_path, timeout=timing)


def test_weirdbasis(tmp_path, timing=30):
    """Runs a runcard with perturbative charm basis but an intrinsic charm theory
    so the run will be stopped by the checks"""
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


@pytest.mark.linux
@pytest.mark.parametrize("runcard", [QUICKNAME_SEQUENTIAL, QUICKNAME_PARALLEL])
def test_multireplica_runs(tmp_path, runcard):
    """`n3fit 1 -r 3`, `n3fit 2 -r 3` and `n3fit 1 -r 2` should give identical results for replica 2"""
    quickcard = f"{runcard}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    options = {'a': '1 -r 3', 'b': '2 -r 3', 'c': '3'}
    jsons = {}
    for name, replicas in options.items():
        path = tmp_path / name
        path.mkdir()
        shutil.copy(quickpath, path)
        sp.run(f"{EXE} {quickcard} {replicas}".split(), cwd=path, check=True)

    for name_1, option_1 in options.items():
        for name_2, option_2 in options.items():
            if name_1 > name_2:
                path_1 = tmp_path / f"{name_1}/{runcard}/nnfit/replica_3/weights.h5"
                path_2 = tmp_path / f"{name_2}/{runcard}/nnfit/replica_3/weights.h5"
                with h5py.File(path_1, 'r') as file_1, h5py.File(path_2, 'r') as file_2:
                    compare_weights(option_1, option_2, file_1, file_2)


def compare_weights(option_1, option_2, file_1, file_2):
    for key in file_1.keys():
        # The bias is initialized to 0 and will have possibly large relative differences
        if 'bias' in key:
            continue
        if isinstance(file_1[key], h5py.Group):
            compare_weights(option_1, option_2, file_1[key], file_2[key])
        else:
            weight_name = file_1[key].name
            err_msg = f"Difference between runs `n3fit {option_1}` and `n3fit {option_2}` in weights {weight_name}"
            np.testing.assert_allclose(
                file_1[key][:], file_2[key][:], rtol=1e-5, atol=1e-5, err_msg=err_msg
            )
