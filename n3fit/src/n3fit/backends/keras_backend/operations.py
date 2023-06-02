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
from tensorflow.keras.layers import subtract as keras_subtract

from tensorflow.keras.layers import Input
from tensorflow.keras import backend as K

from validphys.convolution import OP


def evaluate(tensor):
    """ Evaluate input tensor using the backend """
    return K.eval(tensor)


def as_layer(operation, op_args=None, op_kwargs=None, **kwargs):
    """Wrap any operation as a keras layer

    Note that a layer call argument takes only one argument, therefore
    all extra arguments defining the operation must be given as part
    of `op_args` (a list) and `op_kwargs` (a dict) and will be compiled
    together with the operation

    Parameters
    ----------
        operation: function
            opertion to compute (its first argument must be for a tensor)
        op_args: list
            list of positional arguments for the operation
        op_kwargs: dict
            dict of optional arguments for the operation

    Returns
    -------
        op_layer: layer
            a keras layer that applies the operation upon call
    """
    if op_args is None:
        op_args = []
    if op_kwargs is None:
        op_kwargs = {}

    def apply_op(x):
        return operation(x, *op_args, **op_kwargs)

    op_layer = keras_Lambda(apply_op, **kwargs)
    return op_layer


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

    @tf.function
    def operate_on_tensors(tensor_list):
        return operation(*tensor_list)

    return operate_on_tensors


# f(x: numpy) -> y: tensor
def numpy_to_tensor(ival, **kwargs):
    """
    Make the input into a tensor
    """
    return K.constant(ival, **kwargs)


# f(x: tensor) -> y: tensor
def batchit(x, batch_dimension=0, **kwarg):
    """ Add a batch dimension to tensor x """
    return tf.expand_dims(x, batch_dimension, **kwarg)


# layer generation
def numpy_to_input(
        numpy_array: np.ndarray,
        no_reshape: bool = False,
        name: str = None,
        custom_shape: tuple = None,
        ):
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
        custom_shape: tuple
            To specify a more general shape with None values
    """
    if no_reshape:
        batched_array = numpy_array
        batch_size = numpy_array.shape[0]
        shape = numpy_array.shape[1:]
    else:
        batched_array = np.expand_dims(numpy_array, 0)
        batch_size = 1
        shape = numpy_array.shape
    if custom_shape is not None:
        shape = custom_shape
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

    layer_op = as_layer(lambda inputs: inputs[0] * inputs[1])
    return layer_op(o_list)


def op_gather_keep_dims(tensor, indices, axis=0, **kwargs):
    """A convoluted way of providing ``x[:, indices, :]``

    From TF 2.4 onwards tensorflow is able to understand the syntax above for
    both eager and non-eager tensors
    """
    if indices == -1:
        indices = tensor.shape[axis]-1

    def tmp(x):
        y = tf.gather(x, indices, axis=axis, **kwargs)
        return tf.expand_dims(y, axis=axis)

    layer_op = as_layer(tmp)
    return layer_op(tensor)


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


def stack(tensor_list, axis=0, **kwargs):
    """ Stack a list of tensors
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/stack>`_
    """
    return tf.stack(tensor_list, axis=axis, **kwargs)


def concatenate(tensor_list, axis=-1, target_shape=None, name=None):
    """
    Concatenates a list of numbers or tensor into a bigger tensor
    If the target shape is given, the output is reshaped to said shape
    """
    concatenated_tensor = tf.concat(tensor_list, axis, name=name)
    if target_shape:
        return K.reshape(concatenated_tensor, target_shape)
    else:
        return concatenated_tensor


# Mathematical operations
def pdf_masked_convolution(raw_pdf, basis_mask):
    """Computes a masked convolution of two equal pdfs
    And applies a basis_mask so that only the actually useful values
    of the convolution are returned.

    If training many PDFs at once, it will use as a backend `einsum`, which
    is better suited for running on GPU (but slower on CPU).

    Parameters
    ----------
        pdf: tf.tensor
            rank 4 (batchsize, xgrid, flavours, replicas)
        basis_mask: tf.tensor
            rank  2 tensor (flavours, flavours)
            mask to apply to the pdf convolution

    Return
    ------
        pdf_x_pdf: tf.tensor
            rank3 (len(mask_true), xgrid, xgrid, replicas)
    """
    if raw_pdf.shape[-1] == 1: # only one replica!
        pdf = tf.squeeze(raw_pdf, axis=(0,-1))
        luminosity = tensor_product(pdf, pdf, axes=0)
        lumi_tmp = K.permute_dimensions(luminosity, (3, 1, 2, 0))
        pdf_x_pdf = batchit(boolean_mask(lumi_tmp, basis_mask), -1)
    else:
        pdf = tf.squeeze(raw_pdf, axis=0)  # remove the batchsize
        luminosity = tf.einsum('air,bjr->jibar', pdf, pdf)
        # (xgrid, flavour, xgrid, flavour)
        pdf_x_pdf = boolean_mask(luminosity, basis_mask)
    return pdf_x_pdf


def einsum(equation, *args, **kwargs):
    """
    Computes the tensor product using einsum
    See full `docs <https://www.tensorflow.org/api_docs/python/tf/einsum>`_
    """
    return tf.einsum(equation, *args, **kwargs)


def tensor_product(*args, **kwargs):
    """
    Computes the tensordot product between tensor_x and tensor_y
    See full `docs <https://www.tensorflow.org/api_docs/python/tf/tensordot>`_
    """
    return tf.tensordot(*args, **kwargs)


@tf.function
def pow(tensor, power):
    """
    Computes the power of the tensor
    """
    return tf.pow(tensor, power)


@tf.function(experimental_relax_shapes=True)
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


def split(*args, **kwargs):
    """
    Splits the tensor on the selected axis
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/split>`_
    """
    return tf.split(*args, **kwargs)


def scatter_to_one(values, indices=[[1]], output_dim=14):
    """
    Like scatter_nd initialized to one instead of zero
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/scatter_nd>`_
    """
    ones = np.ones(output_dim, dtype=np.float32)
    return tf.tensor_scatter_nd_update(ones, indices, values)

def scatter_to_zero(values, indices, output_dim):
    """
    Like scatter_nd initialized to zero
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/scatter_nd>`_
    """
    indices = tf.constant([[i] for i in indices])
    updates = tf.constant(values)
    return tf.scatter_nd(indices, updates, [output_dim])


def op_subtract(inputs, **kwargs):
    """
    Computes the difference between two tensors.
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/keras/layers/subtract>`_
    """
    return keras_subtract(inputs, **kwargs)


@tf.function
def backend_function(fun_name, *args, **kwargs):
    """
    Wrapper to call non-explicitly implemented backend functions by name: (``fun_name``)
    see full `docs <https://keras.io/api/utils/backend_utils/>`_ for some possibilities
    """
    fun = getattr(K, fun_name)
    return fun(*args, **kwargs)
