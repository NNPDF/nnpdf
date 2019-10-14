"""
    This module containg a list of useful operations translated in the keras language
    
    The main goal was (is) to implement the NNPDF operations on fktable in the keras
    language (hence the mapping `c_to_py_fun`)

    All operations accept as an input an iterable of keras layers with
    the same output_shape (usually the output of an Observable layer)
    and return a keras operation with, again, the same output shape.
"""

from keras.layers import add as keras_add
from keras.layers import subtract as keras_subtract
from keras.layers import Lambda as keras_Lambda
from keras.layers import multiply as keras_multiply

from keras.layers import Input
from keras import backend as K

import numpy as np

def numpy_to_tensor(ival):
    """
        Make the input into a tensor
    """
    return K.constant(ival)

def numpy_to_input(numpy_array):
    """
        If x is a numpy array, make it into a numpy tensor
        if not just returns the input unchanged
    """
    if isinstance(numpy_array, np.ndarray):
        tensor = K.constant(numpy_array)
        return Input(tensor=tensor)
    else:
        return numpy_array


def c_to_py_fun(op_name, name, default="ADD"):
    """
    Map between the NNPDF operations and the operations defined in this file
    Any new backend must implement such a mapping
    """ # TODO: shouldn't this be outside of the backend folder then
        # surely...
    d = {
        'NULL' : op_null,
        'ADD' : op_add,
        'RATIO' : op_ratio,
        'ASY' : op_asy,
        'SMN' : op_smn,
            }
    if op_name not in d.keys():
        print("Operation name not recognised, defaulting to {0}".format(default))
        return d[default]

    def operation_fun(o_list):
        return d[op_name](o_list, name=name)

    return operation_fun


def op_multiply(o_list, **kwargs):
    """
    Receives a list of layers of the same output size and multiply them element-wise
    """
    return keras_multiply(o_list, **kwargs)


def op_multiply_dim(o_list, **kwargs):
    """
    Bypass in order to multiply two layers with different output dimension
    for instance: (10000 x 14) * (14)
    as the normal keras multiply don't accept it (but somewhow it does accept it doing it like this)
    """
    if len(o_list) != 2:
        raise ValueError(
            "The number of observables is incorrect, operations.py:op_multiply_dim, expected 2, received {0}".format(
                len(o_list)
            )
        )
    create_operation = keras_Lambda(lambda inputs: inputs[0] * inputs[1])
    return create_operation(o_list)


def op_null(o_list, **kwargs):
    """
    Not a compound object, do nothing
    """
    return o_list[0]


def op_add(o_list, **kwargs):
    """
    Sum a list of layers with the same output dim
    """
    return keras_add(o_list, **kwargs)


def op_subtract(o_list, **kwargs):
    """
    Subtract all observables
    """
    return keras_subtract(o_list)


def op_log(o_tensor, **kwargs):
    """
    Computes the logarithm of the input
    """
    return K.log(o_tensor)


def op_ratio(o_list, **kwargs):
    """
    Take the ratio of two observables
    """
    if len(o_list) != 2:
        raise ValueError(
            "The number of observables is incorrect, operations.py:op_ratio, expected 2, received {0}".format(
                len(o_list)
            )
        )

    division_layer = keras_Lambda(lambda inputs: inputs[0] / inputs[1], **kwargs)
    return division_layer(o_list)


def op_asy(o_list, **kwargs):
    """
    Perform the asymmetry operation on two observables
    """
    if len(o_list) != 2:
        raise ValueError(
            "The number of observables is incorrect, operations.py:op_asy, expected 2, received {0}".format(
                len(o_list)
            )
        )

    subtraction = keras_subtract(o_list)
    addition = op_add(o_list)
    return op_ratio([subtraction, addition], **kwargs)


def op_smn(o_list, **kwargs):
    """
    Normalised sum
    """
    if len(o_list) != 4:
        raise ValueError(
            "The number of observables is incorrect, operations.py:op_smn, expected 4, received {0}".format(
                len(o_list)
            )
        )
    numer = op_add(o_list[:2])
    denom = op_add(o_list[2:])
    return op_ratio([numer, denom], **kwargs)
