"""
    This module tests the mathematical functions in the n3fit backend
    and ensures they do the same thing as their numpy counterparts
"""
import operator
import numpy as np
from n3fit.backends import operations as op

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
     - (tensor, array): if passed a tuple (backend tensor, numpy array), uses these
        values as tensor and array inputs for the operations
    """
    if mode == "same":
        tensors = [T1, T2]
        arrays = [ARR1, ARR2]
    elif mode == "diff":
        tensors = [T3, T1]
        arrays = [ARR3, ARR1]
    elif mode == "four":
        tensors = [T1, T2, T1, T1]
        arrays = [ARR1, ARR2, ARR1, ARR1]
    elif mode == "twenty":
        tensors = [T1, T2, T1, T1, T1, T1, T1, T1, T1, T1, T1, T2, T1, T1, T1, T1, T1, T1, T1, T1]
        arrays = [ARR1, ARR2, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR2, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1]
    elif mode == "ten":
        tensors = [T1, T2, T1, T1, T1, T1, T1, T1, T1, T1]
        arrays = [ARR1, ARR2, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1, ARR1]
    elif mode == "single":
        tensors = [T1]
        arrays = [ARR1]
    elif isinstance(mode, tuple):
        tensors = mode[0]
        arrays = mode[1]

    result = backend_op(tensors)
    reference = python_op(*arrays)
    are_equal(result, reference)


# Test NNPDF operations
def test_c_to_py_fun():
    # Null function
    op_null = op.c_to_py_fun("NULL")
    numpy_check(op_null, lambda x: x, "single")
    # Add
    op_add = op.c_to_py_fun("ADD")
    numpy_check(op_add, operator.add)
    # Ratio
    op_rat = op.c_to_py_fun("RATIO")
    numpy_check(op_rat, operator.truediv)
    # ASY
    op_asy = op.c_to_py_fun("ASY")
    reference = lambda x, y: (x - y) / (x + y)
    numpy_check(op_asy, reference)
    # SMN
    op_smn = op.c_to_py_fun("SMN")
    reference = lambda x, y, z, d: (x + y) / (z + d)
    numpy_check(op_smn, reference, "four")
    # COM
    op_com = op.c_to_py_fun("COM")
    reference = lambda x, y, z, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t : (x + y + z + d + e + f + g + h + i + j) / (k + l + m + n + o + p + q + r + s + t)
    numpy_check(op_com, reference, "twenty")
    # SMT
    op_smt = op.c_to_py_fun("SMT")
    reference = lambda x, y, z, d, e, f, g, h, i, j : (x + y + z + d + e + f + g + h + i + j)
    numpy_check(op_smt, reference, "ten")

# Tests operations
def test_op_multiply():
    numpy_check(op.op_multiply, operator.mul)


def test_op_multiply_dim():
    numpy_check(op.op_multiply_dim, operator.mul, mode="diff")


def test_op_log():
    numpy_check(op.op_log, np.log, mode='single')


def test_flatten():
    numpy_check(op.flatten, np.ndarray.flatten, mode=(T3, [ARR3]))


def test_boolean_mask():
    bools = np.random.randint(0, 2, DIM, dtype=bool)
    np_result = ARR1[bools]
    tf_bools = op.numpy_to_tensor(bools)
    tf_result = op.boolean_mask(T1, tf_bools, axis=0)
    are_equal(np_result, tf_result)

def test_tensor_product():
    np_result = np.tensordot(ARR3, ARR1, axes=1)
    tf_result = op.tensor_product(T3, T1, axes=1)
    are_equal(np_result, tf_result)

def test_sum():
    numpy_check(op.sum, np.sum, mode='single')
