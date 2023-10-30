"""
    Test hyperoptimization features
"""
import pathlib
import shutil
import subprocess as sp
import json
import pytest

from numpy.testing import assert_approx_equal
from n3fit.hyper_optimization import rewards


def test_rewards():
    """Ensure that rewards continue doing what they are supposed to do"""
    losses = [0.0, 1.0, 2.0]
    assert_approx_equal(rewards.average(losses), 1.0)
    assert_approx_equal(rewards.best_worst(losses), 2.0)
    assert_approx_equal(rewards.std(losses), 0.816496580927726)


# HYPEROPT_FOLDER = pathlib.Path(__file__).resolve().parent / "hyperopt"
HYPEROPT_FOLDER = pathlib.Path(__file__).with_name("hyperopt")
QUICKNAME = "quickcard"
EXE = "n3fit"
REPLICA = "1"


def load_data(info_file):
    """Loads the information of the fit from the json files"""
    with open(info_file, "r", encoding='utf-8') as file:
        return json.load(file)


def test_restart_from_pickle():
    """Ensure that our hyperopt restart works as expected"""
    # Prepare the run
    quickcard = f"hyper-{QUICKNAME}.yml"
    quickpath = HYPEROPT_FOLDER / quickcard
    # Set up some options
    n_trials_stop = 2
    n_trials_total = 4
    output_restart = (
        HYPEROPT_FOLDER / f"run_{n_trials_stop}_trials_and_then_{n_trials_total}_trials"
    )
    output_direct = HYPEROPT_FOLDER / f"run_{n_trials_total}_trials"

    # cp runcard to tmp folder
    # shutil.copy(quickpath, tmp_path)
    # run some trials for the first time
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_stop} " f"-o {output_restart}".split(),
        check=True,
    )
    # restart and calculate more trials
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_total} "
        f"-o {output_restart} --continue".split(),
        check=True,
    )
    # start again and calculate all trials at once
    sp.run(
        f"{EXE} {quickpath} {REPLICA} --hyperopt {n_trials_total} " f"-o {output_direct}".split(),
        check=True,
    )

    # read up generated json files
    restart_json_path = f"{output_restart}/nnfit/replica_{REPLICA}/tries.json"
    restart_json = load_data(restart_json_path)
    direct_json_path = f"{output_direct}/nnfit/replica_{REPLICA}/tries.json"
    direct_json = load_data(direct_json_path)

    # Remove the 'output_restart' and 'output_direct' directories
    shutil.rmtree(output_restart, ignore_errors=True)
    shutil.rmtree(output_direct, ignore_errors=True)

    # minimum check: the generated list of nested dictionaries have same lenght
    assert len(restart_json) == len(direct_json)

    for i in range(n_trials_total):
        # check that the files share exactly the same hyperopt history
        assert restart_json[i]['misc'] == direct_json[i]['misc']
        assert restart_json[i]['state'] == direct_json[i]['state']
        assert restart_json[i]['tid'] == direct_json[i]['tid']
        assert restart_json[i]['result'] == direct_json[i]['result']
