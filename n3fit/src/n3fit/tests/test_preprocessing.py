import numpy as np

from n3fit.backends import Input, Lambda, MetaModel
from n3fit.backends import operations as op
from n3fit.layers import Preprocessing


def setup_layer(replica_seeds):
    """Setup a layer for testing"""
    # taken from basic runcard
    flav_info = [
        {'fl': 'sng', 'smallx': [1.05, 1.19], 'largex': [1.47, 2.7], 'trainable': False},
        {'fl': 'g', 'smallx': [0.94, 1.25], 'largex': [0.11, 5.87], 'trainable': False},
        {'fl': 'v', 'smallx': [0.54, 0.75], 'largex': [1.15, 2.76], 'trainable': False},
        {'fl': 'v3', 'smallx': [0.21, 0.57], 'largex': [1.35, 3.08]},
        {'fl': 'v8', 'smallx': [0.52, 0.76], 'largex': [0.77, 3.56]},
        {'fl': 't3', 'smallx': [-0.37, 1.52], 'largex': [1.74, 3.39]},
        {'fl': 't8', 'smallx': [0.56, 1.29], 'largex': [1.45, 3.03]},
        {'fl': 'cp', 'smallx': [0.12, 1.19], 'largex': [1.83, 6.7]},
    ]
    prepro = Preprocessing(flav_info=flav_info, replica_seeds=replica_seeds)
    np.random.seed(42)
    test_x = np.random.uniform(size=(1, 4, 1))
    return prepro, test_x


def test_preprocessing():
    """Regression test"""
    prepro, test_x = setup_layer(replica_seeds=[1])
    test_prefactors = [
        [
            [
                [
                    3.7446213e-01,
                    1.9785003e-01,
                    2.7931085e-01,
                    2.0784079e-01,
                    4.5369801e-01,
                    2.7796263e-01,
                    5.4610312e-01,
                    2.4907256e-02,
                ],
                [
                    6.2252983e-04,
                    3.0504008e-05,
                    4.5713778e-03,
                    1.0905267e-03,
                    4.0506415e-02,
                    5.9004971e-05,
                    4.5114113e-03,
                    2.6757403e-09,
                ],
                [
                    4.1631009e-02,
                    1.0586979e-02,
                    8.3202787e-02,
                    4.3506064e-02,
                    2.2559988e-01,
                    1.5161950e-02,
                    1.0105091e-01,
                    1.4808348e-04,
                ],
                [
                    1.1616933e-01,
                    4.2717375e-02,
                    1.5620175e-01,
                    9.7478621e-02,
                    3.2600221e-01,
                    5.8901049e-02,
                    2.1937098e-01,
                    1.8343410e-03,
                ],
            ]
        ]
    ]
    prefactors = prepro(test_x)
    np.testing.assert_allclose(test_prefactors, prefactors, rtol=1e-6)


def test_constraint():
    """Test the constraint"""
    prepro, test_x = setup_layer(replica_seeds=[1, 5])
    prefactors = prepro(test_x)

    # create model that we can train
    x = Input(shape=test_x.shape[1:])
    prefactors = prepro(x)
    scalar = Lambda(lambda x: op.sum(x, axis=(1, 2, 3)))(prefactors)
    model = MetaModel(input_tensors={'x': x}, output_tensors=scalar)
    model.compile(loss='mse', learning_rate=1e-15)

    # Simulate training where weights of replica 1 are updated to violate the constraint
    prepro.weights[0].assign(10.0 * prepro.weights[0])
    weights_before = [w.numpy() for w in prepro.weights]
    # Check that indeed we violate the constraint now
    assert np.any(prepro.weights[0].numpy() > prepro.weights[0].constraint.max_value)

    # Train for one step to let the constraint kick in
    model.fit(test_x, np.array([0.0]), epochs=1)
    weights_after = [w.numpy() for w in prepro.weights]

    # Check that now everything satisfies the constraint again
    for w in prepro.weights:
        if w.trainable:
            assert np.alltrue(w.constraint.min_value <= w.numpy())
            assert np.alltrue(w.numpy() <= w.constraint.max_value)

    # Check that other replicas were not affected
    for wa, wb in zip(weights_after[1:], weights_before[1:]):
        np.testing.assert_allclose(wa, wb, atol=1e-6)
