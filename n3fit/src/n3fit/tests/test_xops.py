"""
    Test the x operations
"""
import numpy as np

from n3fit.layers import xDivide


def test_xdivide_default():
    """Check that the default xDivide works as expected"""
    x_div = xDivide()
    test_input = np.array([1, 2, 3], dtype=np.float32).reshape((1, 3, 1))
    test_output = x_div(test_input)

    expected_output = np.ones(shape=(1, 3, 14))
    default_indices = [3, 4, 5, 6]
    for i in default_indices:
        expected_output[:, :, i] = 1 / test_input[:, :, 0]

    np.testing.assert_allclose(test_output, expected_output, rtol=1e-05)


def test_xdivide_indices():
    """Check that the default xDivide works as expected"""
    custom_indices = [0, 1, 7]
    x_div = xDivide(div_list=custom_indices)
    test_input = np.array([1, 2, 3], dtype=np.float32).reshape((1, 3, 1))
    test_output = x_div(test_input)

    expected_output = np.ones(shape=(1, 3, 14))
    for i in custom_indices:
        expected_output[:, :, i] = 1 / test_input[:, :, 0]

    np.testing.assert_allclose(test_output, expected_output, rtol=1e-05)
