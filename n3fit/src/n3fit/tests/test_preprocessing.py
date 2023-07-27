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
    test_prefactors = np.array(
        [
            [
                [
                    4.65414822e-01,
                    3.86829615e-01,
                    3.37056696e-01,
                    2.07052663e-01,
                    2.88195640e-01,
                    1.44142181e-01,
                    2.54527658e-01,
                    5.00672087e-02,
                ],
                [
                    4.95845545e-03,
                    8.82980123e-04,
                    4.71688760e-03,
                    1.34412316e-03,
                    2.12713517e-03,
                    3.22611211e-03,
                    1.32401168e-04,
                    4.70521542e-08,
                ],
                [
                    9.99829322e-02,
                    4.81918976e-02,
                    8.90727192e-02,
                    4.71593030e-02,
                    6.22668676e-02,
                    5.96290566e-02,
                    2.02597082e-02,
                    5.59301290e-04,
                ],
                [
                    2.06480861e-01,
                    1.27716750e-01,
                    1.73152477e-01,
                    1.01946764e-01,
                    1.33761078e-01,
                    1.03049524e-01,
                    6.74647838e-02,
                    4.97586234e-03,
                ],
            ]
        ]
    )
    prefactors = prepro(test_x)
    assert np.allclose(test_prefactors, prefactors)
