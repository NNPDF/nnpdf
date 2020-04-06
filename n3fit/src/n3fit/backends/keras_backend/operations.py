"""
    This module containg a list of useful operations translated in the keras language

    All operations accept as input an iterable of keras layers or tensors
    and (when necessary) keyword arguments.
    The return operation is always a keras layer (or tensor)

    This includes an implementation of the NNPDF operations on fktable in the keras
    language (hence the mapping `c_to_py_fun`)
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Lambda as keras_Lambda
from tensorflow.keras.layers import multiply as keras_multiply
from tensorflow.keras.layers import Concatenate as keras_concatenate

from tensorflow.keras.layers import Input, Layer
from tensorflow.keras import backend as K

from validphys.convolution import OP


def numpy_to_tensor(ival):
    """
        Make the input into a tensor
    """
    return K.constant(ival)


def batchit(x, batch_dimension=0):
    """ Add a batch dimension to tensor x """
    return tf.expand_dims(x, batch_dimension)

def concatenate_split(splitting_sizes, axis=1):
    """ Generate a pair of concatention and splitting layer
    so that they invert each other

    Parameters
    ----------
        splitting_sizes: list(int)
            size of the output of the split
        axis: int
            axis in which to apply the operation
    """
    concatenation_layer = keras_concatenate(axis=axis)
    splitting_layer = keras_Lambda( lambda x: tf.split(x, splitting_sizes, axis=axis) )
    return concatenation_layer, splitting_layer


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


def c_to_py_fun(op_name, name = "dataset"):
    """
    Map the NNPDF operations to Keras layers
    NNPDF operations are defined in :py:func:`validphys.convolution.OP

    Parameters
    ----------
        op_name: str
            A string defining the operation name
    """
    try:
        operation = OP[op_name]
    except KeyError as e:
        raise ValueError(f"Operation {op_name} not recognised") from e

    # Convert the operation into a lambda layer
    operation_layer = keras_Lambda(lambda x: operation(*x), name=f"op_{name}_{op_name}")
    return operation_layer

def op_subtract(o_list): #TODO to be removed, not used once other PRs are merged
    from tensorflow.keras.layers import subtract
    return subtract(o_list)


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


def op_log(o_tensor, **kwargs):
    """
    Computes the logarithm of the input
    """
    return K.log(o_tensor)
