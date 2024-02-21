"""
    Test hyperoptimization features
"""
import json
import pathlib
import shutil
import subprocess as sp
import time

from numpy.testing import assert_approx_equal

from n3fit.hyper_optimization import rewards


def test_rewards():
    """Ensure that rewards continue doing what they are supposed to do"""
    losses = [0.0, 1.0, 2.0]
    assert_approx_equal(rewards.average(losses), 1.0)
    assert_approx_equal(rewards.best_worst(losses), 2.0)
    assert_approx_equal(rewards.std(losses), 0.816496580927726)


REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
QUICKNAME = "quickcard"
EXE = "n3fit"
REPLICA = "1"


def load_data(info_file):
    """Loads the information of the fit from the json files"""
    with open(info_file, "r", encoding='utf-8') as file:
        return json.load(file)


def test_restart_from_pickle(tmp_path):
    """Ensure that our hyperopt restart works as expected"""
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = REGRESSION_FOLDER / quickcard
    # Set up some options
    n_trials_stop = 2
    n_trials_total = 4
    output_restart = tmp_path / f"run_{n_trials_stop}_trials_and_then_{n_trials_total}_trials"
    output_direct = tmp_path / f"run_{n_trials_total}_trials"

    # cp runcard to tmp folder
    shutil.copy(quickpath, tmp_path)
    # run some trials for the first time
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_stop} " f"-o {output_restart}".split(),
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

    # minimum check: the generated list of nested dictionaries have same lenght
    assert len(restart_json) == len(direct_json)

    for i in range(n_trials_total):
        # check that the files share exactly the same hyperopt history
        assert restart_json[i]['misc'] == direct_json[i]['misc']
        assert restart_json[i]['state'] == direct_json[i]['state']
        assert restart_json[i]['tid'] == direct_json[i]['tid']
        assert restart_json[i]['result'] == direct_json[i]['result']


def start_mongo_database(tmp_path):
    """Creates MongoDB database and returns the Popen object."""
    db_command = ["mongod", "--dbpath", f"{tmp_path}/hyperopt-db"]
    directory_path = f"{tmp_path}/hyperopt-db"
    try:
        # create database directory
        sp.run(["mkdir", "-p", directory_path], check=True)
        # launch database
        process = sp.Popen(db_command, cwd=tmp_path)
        return process
    except (sp.CalledProcessError, OSError) as err:
        msg = f"Error creating directory or executing {db_command}: {err}"
        raise EnvironmentError(msg) from err


def stop_mongod_command(process):
    """Stops the MongoDB database."""
    # directory_path = f"{tmp_path}/hyperopt"
    try:
        # stop mongod command
        process.terminate()
        process.wait()
        # remove database files
        # sp.run(f"rm -r {directory_path} && rm -r {tmp_path}/65*", check=True)
    except (sp.CalledProcessError, OSError) as err:
        msg = f"Error stopping the MongoDB process or removing database files: {err}"
        raise EnvironmentError(msg) from err


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

    # Generate on-the-fly a real MongoDB database
    process = start_mongo_database(tmp_path)

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

    # Stop mongod command
    stop_mongod_command(process)

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
