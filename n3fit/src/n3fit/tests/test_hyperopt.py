"""
    Test hyperoptimization features
"""
import json
import pathlib
import shutil
import subprocess as sp
import tarfile
import time

import numpy as np
from numpy.testing import assert_approx_equal
import pytest

from n3fit.hyper_optimization.rewards import HyperLoss
from n3fit.model_gen import generate_pdf_model
from validphys.loader import Loader


def generate_pdf(seed, num_replicas):
    """Generate generic pdf model."""
    fake_fl = [
        {"fl": i, "largex": [0, 1], "smallx": [1, 2]}
        for i in ["u", "ubar", "d", "dbar", "c", "g", "s", "sbar"]
    ]
    pdf_model = generate_pdf_model(
        nodes=[8],
        activations=["linear"],
        seed=seed,
        num_replicas=num_replicas,
        flav_info=fake_fl,
        fitbasis="FLAVOUR",
    )
    return pdf_model


def get_experimental_data(dataset_name="NMC_NC_NOTFIXED_P_EM-SIGMARED", theoryid=399):
    """Get experimental data set using validphys.

    Returns a list defined by the data set as
    `validphys.core.DataGroupSpec`.
    """
    loader = Loader()
    ds = loader.check_dataset(dataset_name, theoryid=theoryid, cuts="internal", variant="legacy")
    return loader.check_experiment("My DataGroupSpec", [ds])


@pytest.mark.parametrize(
    "loss_type, replica_statistic, expected_per_fold_loss",
    [
        ("chi2", "average", 0.15),
        ("chi2", "best_worst", 0.2),
        ("chi2", "std", 0.05),
        ("phi2", None, None),
    ],
)
def test_compute_per_fold_loss(loss_type, replica_statistic, expected_per_fold_loss):
    """Check that the losses per fold are calculated correctly.

    This example assumes a 2 replica calculation with 3 added penalties.
    """
    # generate 2 replica pdf model
    pdf_model = generate_pdf(seed=[1, 2], num_replicas=2)
    # add 3 penalties for a 2 replica model
    penalties = {
        'saturation': np.array([0.0, 0.0]),
        'patience': np.array([0.0, 0.0]),
        'integrability': np.array([0.0, 0.0]),
    }
    # experimental losses for each replica
    experimental_loss = np.array([0.1, 0.2])
    # get experimental data to compare with
    experimental_data = [get_experimental_data()]

    loss = HyperLoss(loss_type=loss_type, replica_statistic=replica_statistic)

    # calculate statistic loss for one specific fold
    predicted_per_fold_loss = loss.compute_loss(
        penalties, experimental_loss, pdf_model, experimental_data
    )

    # Assert
    if expected_per_fold_loss is not None:
        assert_approx_equal(predicted_per_fold_loss, expected_per_fold_loss)
    else:
        assert predicted_per_fold_loss > 0  # Test for non-negativity
        assert predicted_per_fold_loss.dtype == np.float64  # Test its type
        # Add more property-based tests specific to "phi2" if possible


def test_loss_reduce_over_folds():
    """Ensure that the hyper loss statistics over all folds are calculated correctly."""
    # define losses for 3 folds
    losses = np.array([1.0, 2.0, 3.0])

    loss_average = HyperLoss(fold_statistic="average")
    assert_approx_equal(loss_average.reduce_over_folds(losses), 2.0)

    loss_best_worst_best_worst = HyperLoss(fold_statistic="best_worst")
    assert_approx_equal(loss_best_worst_best_worst.reduce_over_folds(losses), 3.0)

    loss_std = HyperLoss(fold_statistic="std")
    assert_approx_equal(loss_std.reduce_over_folds(losses), 0.816496580927726)


REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
EXE = "n3fit"
REPLICA = "1"


def load_data(info_file):
    """Loads the information of the fit from the json files"""
    with open(info_file, "r", encoding='utf-8') as file:
        return json.load(file)


def test_restart_from_pickle(tmp_path):
    """Ensure that after a hyperopt restart, the testing continues
    from the same point.
    The test is set up so that it does one trial, then stops, then a second one
    And then this is compared with two trials one after the other.

    The test checks that the starting point of the second trial is the same in both cases
    """
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard

    # Set the test up so that it does one trial, then stops, then does another one
    # and then we do two
    n_trials_stop = 1
    n_trials_total = 2
    output_restart = tmp_path / f"run_{n_trials_stop}_trials_and_then_{n_trials_total}_trials"
    output_direct = tmp_path / f"run_{n_trials_total}_trials"

    # cp runcard to tmp folder
    shutil.copy(quickpath, tmp_path)
    # run some trials for the first time
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_stop} -o {output_restart}".split(),
        cwd=tmp_path,
        check=True,
    )
    # restart and calculate more trials
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_total} "
        f"-o {output_restart} --restart".split(),
        cwd=tmp_path,
        check=True,
    )
    # start again and calculate all trials at once
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_total} " f"-o {output_direct}".split(),
        cwd=tmp_path,
        check=True,
    )

    # read up generated json files
    restart_json_path = f"{output_restart}/nnfit/replica_{REPLICA}/tries.json"
    restart_json = load_data(restart_json_path)
    direct_json_path = f"{output_direct}/nnfit/replica_{REPLICA}/tries.json"
    direct_json = load_data(direct_json_path)

    # minimum check: the generated list of nested dictionaries have same length
    assert len(restart_json) == len(direct_json)

    for i in range(n_trials_total):
        # check that the files share exactly the same hyperopt history
        assert restart_json[i]['misc'] == direct_json[i]['misc']
        assert restart_json[i]['state'] == direct_json[i]['state']
        assert restart_json[i]['tid'] == direct_json[i]['tid']
        assert restart_json[i]['misc']['idxs'] == direct_json[i]['misc']['idxs']
    # Note that it doesn't check the final loss of the second trial

@pytest.mark.linux
def test_parallel_hyperopt(tmp_path):
    """Ensure that the parallel implementation of hyperopt with MongoDB works as expected."""
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard

    # Define number of trials and number of mongo-workers to launch
    n_trials = 6
    n_mongo_workers = 3

    # Set up output directories
    output_sequential = tmp_path / "run_hyperopt_sequential"
    output_parallel = tmp_path / "run_hyperopt_parallel"

    # cp runcard to tmp folder
    shutil.copy(quickpath, tmp_path)

    # Run hyperopt sequentially
    start_time = time.time()
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials} " f"-o {output_sequential}".split(),
        cwd=tmp_path,
        check=True,
    )
    end_time = time.time()
    sequential_run_time = end_time - start_time

    # Run hyperopt in parallel
    start_time = time.time()
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials} "
        f"--parallel-hyperopt --num-mongo-workers {n_mongo_workers} "
        f"-o {output_parallel}".split(),
        cwd=tmp_path,
        check=True,
    )
    end_time = time.time()
    parallel_run_time = end_time - start_time

    # Read up generated json files
    sequential_json_path = f"{output_sequential}/nnfit/replica_{REPLICA}/tries.json"
    sequential_json = load_data(sequential_json_path)
    parallel_json_path = f"{output_parallel}/nnfit/replica_{REPLICA}/tries.json"
    parallel_json = load_data(parallel_json_path)

    # Check that the parallel run time is lower than the sequential one
    assert parallel_run_time < sequential_run_time

    # Check that the final json files have the same number of trials
    assert len(parallel_json) == len(sequential_json)

    for i in range(n_trials):
        # Check that the files share the same content
        assert len(parallel_json[i]['misc']) == len(sequential_json[i]['misc'])
        assert len(parallel_json[i]['result']) == len(sequential_json[i]['result'])
        # Note: cannot check that they share exactly the same history
        # as the hyperopt algorithm depends on the results from previous runs
        # which is obviously different between parallel and sequential runs


def clean_up_database(tmp_path, database_name):
    """Stops the MongoDB database."""
    directory_path = f"{tmp_path}/{database_name}"
    try:
        sp.run(f"rm -r {directory_path}", shell=True, check=True)
    except (sp.CalledProcessError, OSError) as err:
        msg = f"Error cleaning up database: {err}"
        raise EnvironmentError(msg) from err


def get_tar_size(filetar):
    """Returns the size of a tar file."""

    def tar_size(tar):
        return sum(member.size for member in tar.getmembers())

    with tarfile.open(filetar, 'r') as tar:
        size = tar_size(tar)
    return size


def test_restart_from_tar(tmp_path):
    """Ensure that our parallel hyperopt restart works as expected."""
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard

    # Set up some options
    n_mongo_workers = 3
    n_trials_stop = 3
    n_trials_total = 6
    output = tmp_path / "output"
    database_name = f"hyperopt-db-{output.name}"

    # cp runcard to tmp folder
    shutil.copy(quickpath, tmp_path)
    # run some trials for the first time
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_stop} "
        f"--parallel-hyperopt --num-mongo-workers {n_mongo_workers} "
        f"-o {output}".split(),
        cwd=tmp_path,
        check=True,
    )
    json_path = f"{output}/nnfit/replica_{REPLICA}/tries.json"
    tar_name = f"{output}/nnfit/replica_{REPLICA}/{database_name}.tar.gz"
    initial_json = load_data(json_path)
    initial_tar_size = get_tar_size(tar_name)

    # just in case, remove old database files to ensure that the restart occurs via tar file
    clean_up_database(tmp_path, database_name)

    # restart and calculate more trials
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_total} "
        f"--parallel-hyperopt --num-mongo-workers {n_mongo_workers} "
        f"-o {output} --restart".split(),
        cwd=tmp_path,
        check=True,
    )
    final_json = load_data(json_path)
    final_tar_size = get_tar_size(tar_name)

    # check if the calculations went well
    assert len(initial_json) == n_trials_stop
    assert len(final_json) == n_trials_total

    # check if the tar files were generated correctly
    assert tarfile.is_tarfile(tar_name) is True

    # check if the final tar file was updated after restart
    assert final_tar_size > initial_tar_size

    for i in range(n_trials_stop):
        # check that the json files share exactly the same hyperopt history until the restart
        assert initial_json[i]['misc'] == final_json[i]['misc']
        assert initial_json[i]['state'] == final_json[i]['state']
        assert initial_json[i]['tid'] == final_json[i]['tid']
        assert initial_json[i]['result'] == final_json[i]['result']
