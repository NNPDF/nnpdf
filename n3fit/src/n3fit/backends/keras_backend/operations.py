"""
    This module containg a list of useful operations translated in the keras language

    All operations accept as input an iterable of keras layers or tensors
    and (when necessary) keyword arguments.
    The return operation is always a keras layer (or tensor)

    This includes an implementation of the NNPDF operations on fktable in the keras
    language (hence the mapping `c_to_py_fun`)
"""

import tensorflow as tf
from tensorflow.keras.layers import add as keras_add
from tensorflow.keras.layers import subtract as keras_subtract
from tensorflow.keras.layers import Lambda as keras_Lambda
from tensorflow.keras.layers import multiply as keras_multiply

from tensorflow.keras.layers import Input, Layer
from tensorflow.keras import backend as K

import numpy as np


def numpy_to_tensor(ival):
    """
        Make the input into a tensor
    """
    return K.constant(ival)


def batchit(x):
    """ Add a batch dimension to tensor x """
    return tf.expand_dims(x, 0)


def numpy_to_input(numpy_array, no_reshape=False):
    """
    Takes a numpy array and generates a Input layer.
    By default it adds a batch dimension (of size 1) so that the shape of the layer
    is that of the array

    Parameters
    ----------
        numpy_array: np.ndarray
        no_reshape: bool
            if true, don't add batch dimension, take the first dimension of the array as the batch
    """
    if isinstance(numpy_array, np.ndarray):
        if no_reshape:
            batched_array = numpy_array
            batch_size = numpy_array.shape[0]
            shape = numpy_array.shape[1:]
        else:
            batched_array = np.expand_dims(numpy_array, 0)
            batch_size = 1
            shape = numpy_array.shape
        input_layer = Input(batch_size=batch_size, shape=shape)
        input_layer.tensor_content = batched_array
        input_layer.original_shape = no_reshape
        return input_layer
    else:
        return numpy_array


def evaluate(tensor):
    """ Evaluate input tensor using the backend """
    return K.eval(tensor)


def c_to_py_fun(op_name, name, default="ADD"):
    """
    Map between the NNPDF operations and the operations defined in this file
    Any new backend must implement such a mapping
    """  # TODO: shouldn't this be outside of the backend folder then
    # surely...
    d = {
        "NULL": op_null,
        "ADD": op_add,
        "RATIO": op_ratio,
        "ASY": op_asy,
        "SMN": op_smn,
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

def m_tensor_ones_like(*args, **kwargs):
    return K.ones_like(*args, **kwargs)

