"""
    This module contains the list of operations that can be used within the
    ``call`` method of the ``n3fit`` layers as well as operations that can
    act on layers.

    This includes an implementation of the NNPDF operations on fktable in the keras
    language (with the mapping ``c_to_py_fun``) into Keras ``Lambda`` layers.

    Tensor operations are compiled through the @tf.function decorator for optimization

    The rest of the operations in this module are divided into four categories:
    numpy to tensor:
        Operations that take a numpy array and return a tensorflow tensor
    layer to layer:
        Operations that take a layer and return another layer
    tensor to tensor:
        Operations that take a tensor and return a tensor
    layer generation:
        Instanciate a layer to be applied by the calling function

    Some of these are just aliases to the backend (tensorflow or Keras) operations
    Note that tensor operations can also be applied to layers as the output of a layer is a tensor
    equally operations are automatically converted to layers when used as such.
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Lambda as keras_Lambda
from tensorflow.keras.layers import multiply as keras_multiply
from tensorflow.keras.layers import Concatenate as keras_concatenate

from tensorflow.keras.layers import Input
from tensorflow.keras import backend as K

from validphys.convolution import OP


def evaluate(tensor):
    """ Evaluate input tensor using the backend """
    return K.eval(tensor)


# NNPDF operations
def c_to_py_fun(op_name, name="dataset"):
    """
    Map the NNPDF operations to Keras layers
    NNPDF operations are defined in :py:func:`validphys.convolution.OP`

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


# f(x: numpy) -> y: tensor
def numpy_to_tensor(ival, **kwargs):
    """
        Make the input into a tensor
    """
    return K.constant(ival, **kwargs)


# f(x: tensor) -> y: tensor
@tf.function
def batchit(x, batch_dimension=0):
    """ Add a batch dimension to tensor x """
    return tf.expand_dims(x, batch_dimension)


# layer generation
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
    splitting_layer = keras_Lambda(lambda x: tf.split(x, splitting_sizes, axis=axis))
    return concatenation_layer, splitting_layer


# layer generation
def numpy_to_input(numpy_array, no_reshape=False, name=None):
    """
    Takes a numpy array and generates a Input layer.
    By default it adds a batch dimension (of size 1) so that the shape of the layer
    is that of the array

    Parameters
    ----------
        numpy_array: np.ndarray
        no_reshape: bool
            if true, don't add batch dimension, take the first dimension of the array as the batch
        name: bool
            name to give to the layer
    """
    if no_reshape:
        batched_array = numpy_array
        batch_size = numpy_array.shape[0]
        shape = numpy_array.shape[1:]
    else:
        batched_array = np.expand_dims(numpy_array, 0)
        batch_size = 1
        shape = numpy_array.shape
    input_layer = Input(batch_size=batch_size, shape=shape, name=name)
    input_layer.tensor_content = batched_array
    input_layer.original_shape = no_reshape
    return input_layer


#
# Layer to Layer operations
#
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


#
# Tensor operations
# f(x: tensor[s]) -> y: tensor
#

# Generation operations
# generate tensors of given shape/content
@tf.function
def tensor_ones_like(*args, **kwargs):
    """
    Generates a tensor of ones of the same shape as the input tensor
    See full `docs <https://www.tensorflow.org/api_docs/python/tf/keras/backend/ones_like>`_
    """
    return K.ones_like(*args, **kwargs)


@tf.function
def many_replication(grid, replications, axis=0, **kwargs):
    """
    Generates a tensor with one extra dimension:
        a repetition of "grid" n times along the given axis
    from keras documentation:
    If x has shape (s1, s2, s3) and axis is 1, the output will have shape (s1, s2 * rep, s3)
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/keras/backend/repeat_elements>`_
    """
    return K.repeat_elements(grid, rep=replications, axis=axis, **kwargs)


# Property operations
# modify properties of the tensor like the shape or elements it has
@tf.function
def flatten(x):
    """ Flatten tensor x """
    return tf.reshape(x, (-1,))


def boolean_mask(*args, **kwargs):
    """
    Applies a boolean mask to a tensor

    Relevant parameters: (tensor, mask, axis=None)
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/boolean_mask>`_.
    """
    return tf.boolean_mask(*args, **kwargs)


@tf.function
def transpose(tensor, **kwargs):
    """
    Transpose a layer,
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/keras/backend/transpose>`_
    """
    return K.transpose(tensor, **kwargs)


def concatenate(tensor_list, axis=-1, target_shape=None):
    """
    Concatenates a list of numbers or tenosr into a bigger tensor
    If the target shape is given, the output is reshaped to said shape
    """
    concatenated_tensor = K.concatenate(tensor_list, axis=axis)
    if target_shape:
        return K.reshape(concatenated_tensor, target_shape)
    else:
        return concatenated_tensor


# Mathematical operations
def pdf_masked_convolution(raw_pdf, basis_mask):
    """ Computes a masked convolution of two equal pdfs
    And applies a basis_mask so that only the actually useful values
    of the convolution are returned

    Parameters
    ----------
        pdf: tf.tensor
            rank 3 (batchsize, xgrid, flavours)
        basis_mask: tf.tensor
            rank  2 tensor (flavours, flavours)
            mask to apply to the pdf convolution

    Return
    ------
        pdf_x_pdf: tf.tensor
            rank3 (len(mask_true), xgrid, xgrid)
    """
    pdf = tf.squeeze(raw_pdf, axis=0)  # remove the batchsize
    luminosity = tensor_product(pdf, pdf, axes=0)
    # (xgrid, flavour, xgrid, flavour)
    # reshape to put the flavour indices at the beginning to apply mask
    lumi_tmp = K.permute_dimensions(luminosity, (3, 1, 2, 0))
    pdf_x_pdf = boolean_mask(lumi_tmp, basis_mask)
    return pdf_x_pdf


def tensor_product(*args, **kwargs):
    """
    Computes the tensordot product between tensor_x and tensor_y
    See full `docs <https://www.tensorflow.org/api_docs/python/tf/tensordot>`_
    """
    return tf.tensordot(*args, **kwargs)


@tf.function
def op_log(o_tensor, **kwargs):
    """
    Computes the logarithm of the input
    """
    return K.log(o_tensor)


@tf.function
def sum(*args, **kwargs):
    """
    Computes the sum of the elements of the tensor
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/keras/backend/sum>`_
    """
    return K.sum(*args, **kwargs)
