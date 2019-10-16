"""
    This module tests the mathematical functions in the n3fit backend
    and ensures they do the same thing as their numpy counterparts
"""
import operator
import numpy as np
from n3fit.backends import operations as op
from n3fit.backends import losses

# General parameters
DIM = 7
THRESHOLD = 1e-6

# Arrays to be used during testing
ARR1 = np.random.rand(DIM)
ARR2 = np.random.rand(DIM)
ARR3 = np.random.rand(DIM + 1, DIM)
C = np.random.rand(DIM, DIM)
INVCOVMAT = np.linalg.inv(C @ C.T)

# Backend-tensors to be used during testing
T1 = op.numpy_to_tensor(ARR1)
T2 = op.numpy_to_tensor(ARR2)
T3 = op.numpy_to_tensor(ARR3)


def are_equal(result, reference, threshold=THRESHOLD):
    """ checks the difference between array `reference` and tensor `result` is
    below `threshold` for all elements """
    res = op.evaluate(result)
    assert np.allclose(res, reference, atol=threshold)


def numpy_check(backend_op, python_op, mode="same"):
    """ Receives a backend operation (`backend_op`) and a python operation
    `python_op` and asserts that, applied to two random arrays, the result
    is the same.
    The option `mode` selects the two arrays to be tested and accepts the following
    options:
     - `same` (default): two arrays of the same dimensionality
     - `diff`: first array has one extra dimension that second array
     - `single`: only one array enters the operation
    """
    if mode == "same":
        tensors = [T1, T2]
        arrays = [ARR1, ARR2]
    elif mode == "diff":
        tensors = [T3, T1]
        arrays = [ARR3, ARR1]
    elif mode == "single":
        tensors = T1
        arrays = [ARR1]
    result = backend_op(tensors)
    reference = python_op(*arrays)
    are_equal(result, reference)


# Tests operations
def test_op_multiply():
    numpy_check(op.op_multiply, operator.mul)


def test_op_multiply_dim():
    numpy_check(op.op_multiply_dim, operator.mul, mode="diff")


def test_op_add():
    numpy_check(op.op_add, operator.add)


def test_op_subtract():
    numpy_check(op.op_subtract, operator.sub)


def test_op_log():
    numpy_check(op.op_log, np.log, mode="single")


def test_op_ratio():
    numpy_check(op.op_ratio, operator.truediv)


def test_op_asy():
    result = op.op_asy([T1, T2])
    reference = (ARR1 - ARR2) / (ARR1 + ARR2)
    are_equal(result, reference)


def test_op_smn():
    result = op.op_smn([T1, T2, T1, T1])
    reference = (ARR1 + ARR2) / (ARR1 + ARR1)
    are_equal(result, reference)


# Tests loss functions
def test_l_invcovmat():
    loss_f = losses.l_invcovmat(INVCOVMAT)
    result = loss_f(T1, T2)
    y = ARR1 - ARR2
    tmp = np.dot(INVCOVMAT, y)
    reference = np.dot(y, tmp)
    are_equal(result, reference)


def test_l_positivity():
    alpha = 1e-7
    loss_f = losses.l_positivity(alpha=alpha)
    result = loss_f(0.0, T1)

    def elu_sum(yarr_in):
        yarr = -yarr_in
        res = 0.0
        for y in yarr:
            if y > 0:
                res += y
            else:
                res += alpha * (np.exp(y) - 1)
        return res

    reference = elu_sum(ARR1)
    are_equal(result, reference)
