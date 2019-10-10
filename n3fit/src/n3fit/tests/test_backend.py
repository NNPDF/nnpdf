"""
    This module tests the mathematical functions in the n3fit backend
    and ensures they do the same thing as their numpy counterparts
"""
import operator
import numpy as np
from keras import backend as K
from n3fit.backends import operations as op

DIM = 7
ARR1 = np.random.rand(DIM)
ARR2 = np.random.rand(DIM)
ARR3 = np.random.rand(DIM+1, DIM)
T1 = op.numpy_to_tensor(ARR1)
T2 = op.numpy_to_tensor(ARR2)
T3 = op.numpy_to_tensor(ARR3)
THRESHOLD = 1e-6

def are_equal(result, reference, threshold = THRESHOLD):
    res = K.eval(result)
    total = res-reference
    assert np.sum(total) < threshold

def numpy_check(backend_op, python_op, mode = 0):
    if mode == 0:
        tensors = [T1, T2]
        arrays = [ARR1, ARR2]
    elif mode == 1:
        tensors = [T3, T1]
        arrays = [ARR3, ARR1]
    elif mode == 2:
        tensors = T1
        arrays = [ARR1]
    result = backend_op(tensors)
    reference = python_op(*arrays)
    are_equal(result, reference)

def test_op_multiply():
    numpy_check(op.op_multiply, operator.mul)

def test_op_multiply_dim():
    numpy_check(op.op_multiply_dim, operator.mul, mode=1)

def test_op_add():
    numpy_check(op.op_add, operator.add)

def test_op_subtract():
    numpy_check(op.op_subtract, operator.sub)

def test_op_log():
    numpy_check(op.op_log, np.log, mode = 2)

def test_op_ratio():
    numpy_check(op.op_ratio, operator.truediv)

def test_op_asy():
    result = op.op_asy([T1,T2])
    reference = (ARR1-ARR2)/(ARR1+ARR2)
    are_equal(result, reference)

def test_op_smn():
    result = op.op_smn([T1,T2,T1,T1])
    reference = (ARR1+ARR2)/(ARR1+ARR1)
    are_equal(result, reference)


if __name__ == "__main__":
    import ipdb
    ipdb.set_trace()

