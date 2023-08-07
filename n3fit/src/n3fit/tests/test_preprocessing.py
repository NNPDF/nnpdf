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
    prepro = Preprocessing(flav_info=flav_info, seed=1)
    np.random.seed(42)
    test_x = np.random.uniform(size=(1, 4, 1))
    test_prefactors = [
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
    prefactors = prepro(test_x)
    np.testing.assert_allclose(test_prefactors, prefactors)
