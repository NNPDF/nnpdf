"""
    Test hyperoptimization features
"""
import os
import random as rn

import numpy as np
from numpy.testing import assert_approx_equal
import pytest
import tensorflow as tf

from n3fit.hyper_optimization.rewards import HyperLoss
from n3fit.model_gen import pdfNN_layer_generator
from validphys.api import API


def set_initial_state():
    """
    This function sets the initial internal state for the different components of n3fit.
    """
    seed = 1
    os.environ.setdefault("HYPEROPT_FMIN_SEED", str(seed))
    np.random.seed(seed)
    use_seed = seed  # np.random.randint(0, pow(2, 31))
    rn.seed(use_seed)
    tf.random.set_seed(use_seed)


def generate_pdf(seeds):
    """Generate generic pdf model."""
    fake_fl = [
        {"fl": i, "largex": [0, 1], "smallx": [1, 2]}
        for i in ["u", "ubar", "d", "dbar", "c", "g", "s", "sbar"]
    ]
    set_initial_state()
    pdf_model = pdfNN_layer_generator(
        nodes=[8], activations=["linear"], seed=seeds, flav_info=fake_fl, fitbasis="FLAVOUR"
    )

    for meta_model in pdf_model:
        for var in meta_model.variables:
            print(f"{var.name}: {var.numpy()}")
    return pdf_model


def get_experimental_data(dataset_name="NMC", theoryid=400):
    """Get experimental data set using validphys API.

    Returns a tuple defined by the data set as
    `validphys.core.DataSetSpec` and associated covariant matrix.
    """
    exp_data_set_spec = API.dataset(
        dataset_input={"dataset": dataset_name}, theoryid=theoryid, use_cuts="internal"
    )
    covmat = API.covariance_matrix(
        dataset_input={"dataset": dataset_name}, theoryid=theoryid, use_cuts="internal"
    )
    return (exp_data_set_spec, covmat)


@pytest.mark.parametrize(
    "loss_type, replica_statistic, expected_per_fold_loss",
    [
        ("chi2", "average", 0.15),
        ("chi2", "best_worst", 0.2),
        ("chi2", "std", 0.05),
        ("phi2", None, 1402.896757483432),
    ],
)
def test_compute_per_fold_loss(loss_type, replica_statistic, expected_per_fold_loss):
    """Check that the losses per fold are calculated correctly.

    This example assumes a 2 replica calculation with 3 added penalties.
    """
    # generate 2 replica pdf model
    pdf_models = generate_pdf(seeds=[0, 1])
    # add 3 penalties for a 2 replica model
    penalties = [np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0])]
    # experimental losses for each replica
    experimental_loss = np.array([0.1, 0.2])
    # get experimental data to compare with
    experimental_data = [get_experimental_data()]

    loss = HyperLoss(loss_type=loss_type, replica_statistic=replica_statistic)

    # calculate statistic loss for one specific fold
    predicted_per_fold_loss = loss.compute_loss(
        penalties, experimental_loss, pdf_models, experimental_data
    )

    assert_approx_equal(predicted_per_fold_loss, expected_per_fold_loss)


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
