from n3fit.layers import Preprocessing
import tensorflow as tf
import numpy as np


def test_preprocessing():
    """Regression test"""
    # taken from basic runcard
    flav_info = [
        {'fl': 'sng', 'smallx': [1.05, 1.19], 'largex': [1.47, 2.7]},  # 'trainable': False},
        {'fl': 'g',   'smallx': [0.94, 1.25], 'largex': [0.11, 5.87]},  # 'trainable': False},
        {'fl': 'v',   'smallx': [0.54, 0.75], 'largex': [1.15, 2.76]},  # 'trainable': False},
        {'fl': 'v3',  'smallx': [0.21, 0.57], 'largex': [1.35, 3.08]},
        {'fl': 'v8',  'smallx': [0.52, 0.76], 'largex': [0.77, 3.56]},
        {'fl': 't3',  'smallx': [-0.37, 1.52], 'largex': [1.74, 3.39]},
        {'fl': 't8',  'smallx': [0.56, 1.29], 'largex': [1.45, 3.03]},
        {'fl': 'cp',  'smallx': [0.12, 1.19], 'largex': [1.83, 6.7]},
    ]
    prepro = Preprocessing(flav_info=flav_info, seed=0)
    test_x = tf.random.uniform(shape=(1, 4, 1), seed=42)
    test_prefactors = tf.constant([[
        [4.28409985e-04, 2.33196355e-02, 2.19706148e-02, 3.12583987e-04, 1.52197943e-04, 6.52114686e-05, 1.36401606e-04, 1.18874144e-04],
        [5.75779974e-02, 2.65057504e-01, 2.13224307e-01, 3.90783511e-02, 3.33565399e-02, 1.79548096e-02, 3.96815017e-02, 2.73224674e-02],
        [1.78151987e-02, 1.46395922e-01, 1.27500281e-01, 1.30333444e-02, 9.50834434e-03, 5.17322961e-03, 1.01000341e-02, 7.87989516e-03],
        [2.80656964e-02, 1.83944702e-01, 1.56266689e-01, 2.01078728e-02, 1.55397197e-02, 8.49795807e-03, 1.71364844e-02, 1.28620844e-02],
                                  ]])
    prefactors = prepro(test_x)
    assert np.allclose(test_prefactors, prefactors)

