import numpy as np

from n3fit.layers import Preprocessing


def test_preprocessing():
    """Regression test"""
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
    prepro = Preprocessing(flav_info=flav_info, seed=0)
    np.random.seed(42)
    test_x = np.random.uniform(size=(1, 4, 1))
    test_prefactors = [
        [
            3.2668063e-01,
            6.7284244e-01,
            3.9915803e-01,
            1.5305418e-01,
            1.8249598e-01,
            7.2065055e-02,
            3.1497714e-01,
            1.3243365e-01,
        ],
        [
            4.6502685e-04,
            2.4272988e-02,
            2.2857290e-02,
            3.3989837e-04,
            1.6686485e-04,
            7.2010938e-05,
            1.4990549e-04,
            1.3058851e-04,
        ],
        [
            3.5664584e-02,
            2.0763232e-01,
            1.7362502e-01,
            2.5178187e-02,
            2.0087767e-02,
            1.0968268e-02,
            2.2659913e-02,
            1.6590673e-02,
        ],
        [
            1.0150098e-01,
            3.5557139e-01,
            2.6867521e-01,
            6.4237528e-02,
            5.9941418e-02,
            3.0898115e-02,
            7.7342905e-02,
            4.8169162e-02,
        ],
    ]
    prefactors = prepro(test_x)
    np.testing.assert_allclose(test_prefactors, prefactors)
