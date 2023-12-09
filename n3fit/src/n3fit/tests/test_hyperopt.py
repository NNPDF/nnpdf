"""
    Test hyperoptimization features
"""
import json
import pathlib
import random as rn
import shutil
import subprocess as sp

import numpy as np
from numpy.testing import assert_approx_equal
import pytest
import tensorflow as tf

from n3fit.backends import clear_backend_state
from n3fit.hyper_optimization.rewards import HyperLoss
from n3fit.model_gen import generate_pdf_model
from validphys.loader import Loader


def set_initial_state(seed=1):
    """
    This function sets the initial internal state for the different components of n3fit.

    Important to warrant that pdf_models are always generated with the same parameters.
    """
    np.random.seed(seed)
    rn.seed(seed)
    clear_backend_state()
    tf.random.set_seed(seed)


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


def get_experimental_data(dataset_name="NMC", theoryid=400):
    """Get experimental data set using validphys.

    Returns a list defined by the data set as
    `validphys.core.DataGroupSpec`.
    """
    loader = Loader()
    ds = loader.check_dataset(dataset_name, theoryid=theoryid, cuts="internal")
    return loader.check_experiment("My DataGroupSpec Name", [ds])


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
    set_initial_state()
    pdf_model = generate_pdf(seed=0, num_replicas=2)
    # add 3 penalties for a 2 replica model
    penalties = [np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0])]
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

    # minimum check: the generated list of nested dictionaries have same length
    assert len(restart_json) == len(direct_json)

    for i in range(n_trials_total):
        # check that the files share exactly the same hyperopt history
        assert restart_json[i]['misc'] == direct_json[i]['misc']
        assert restart_json[i]['state'] == direct_json[i]['state']
        assert restart_json[i]['tid'] == direct_json[i]['tid']
        assert restart_json[i]['result'] == direct_json[i]['result']
