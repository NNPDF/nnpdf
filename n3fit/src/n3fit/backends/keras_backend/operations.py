"""
This module contains the list of operations that can be used within the
``call`` method of the ``n3fit`` layers as well as operations that can
act on layers.

This includes an implementation of the NNPDF operations on fktable in the keras
language (with the mapping ``c_to_py_fun``) into Keras ``Lambda`` layers.

The rest of the operations in this module are divided into four categories:
numpy to tensor:
    Operations that take a numpy array and return a tensorflow tensor
layer to layer:
    Operations that take a layer and return another layer
tensor to tensor:
    Operations that take a tensor and return a tensor
layer generation:
    Instanciate a layer to be applied by the calling function

Most of the operations in this module are just aliases to the backend
(Keras in this case) so that, when implementing new backends, it is clear
which operations may need to be overwritten.
For a few selected operations, a more complicated wrapper to e.g., make
them into layers or apply some default, is included.

Note that tensor operations can also be applied to layers as the output of a layer is a tensor
equally operations are automatically converted to layers when used as such.
"""

from keras import backend as K
from keras import ops as Kops
from keras.layers import ELU, Input
from keras.layers import Lambda as keras_Lambda
import numpy as np

from validphys.convolution import OP

# The following operations are either loaded directly from keras and exposed here
# or the name is change slightly (usually for historical or collision reasons,
# e.g., our logs are always logs or we were using the tf version in the past)

# isort: off
from keras.ops import (
    absolute,
    clip,
    einsum,
    expand_dims,
    leaky_relu,
    reshape,
    repeat,
    split,
    sum,
    tanh,
    transpose,
)
from keras.ops import log as op_log
from keras.ops import power as pow
from keras.ops import take as gather
from keras.ops import tensordot as tensor_product
from keras.layers import multiply as op_multiply
from keras.layers import subtract as op_subtract

# isort: on

# Backend dependent functions and operations
if K.backend() == "torch":
    tensor_to_numpy_or_python = lambda x: x.detach().cpu().numpy()
    decorator_compiler = lambda f: f
elif K.backend() == "jax":
    tensor_to_numpy_or_python = lambda x: np.array(x.block_until_ready())
    decorator_compiler = lambda f: f
elif K.backend() == "tensorflow":
    tensor_to_numpy_or_python = lambda x: x.numpy()
    import tensorflow as tf

    decorator_compiler = tf.function

dict_to_numpy_or_python = lambda ret: {k: tensor_to_numpy_or_python(i) for k, i in ret.items()}
variable_to_numpy = lambda x: x.numpy()


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

    def operate_on_tensors(tensor_list):
        return operation(*tensor_list)

    return operate_on_tensors


# f(x: numpy) -> y: tensor
def numpy_to_tensor(ival, **kwargs):
    """
    Make the input into a tensor
    """
    if (dtype := kwargs.get("dtype", None)) is not bool:
        dtype = K.floatx()
    return Kops.cast(ival, dtype)


# f(x: tensor) -> y: tensor
def batchit(x, batch_dimension=0, **kwarg):
    """Add a batch dimension to tensor x"""
    return Kops.expand_dims(x, batch_dimension, **kwarg)


# layer generation
def numpy_to_input(numpy_array, name=None):
    """
    Takes a numpy array and generates an Input layer with the same shape,
    but with a batch dimension (of size 1) added.

    Parameters
    ----------
        numpy_array: np.ndarray
        name: str
            name to give to the layer
    """
    batched_array = np.expand_dims(numpy_array, 0)
    shape = list(numpy_array.shape)
    # set the number of gridpoints to None, otherwise shapes don't show in model.summary
    shape[0] = None

    input_layer = Input(batch_size=1, shape=shape, name=name)
    input_layer.tensor_content = numpy_to_tensor(batched_array)
    return input_layer


def op_gather_keep_dims(tensor, indices, axis=0, **kwargs):
    """A convoluted way of providing ``x[:, indices, :]``

    From TF 2.4 onwards tensorflow is able to understand the syntax above for
    both eager and non-eager tensors
    """
    if indices == -1:
        indices = tensor.shape[axis] - 1

    def tmp(x):
        y = gather(x, indices, axis=axis)
        return Kops.expand_dims(y, axis=axis)

    layer_op = as_layer(tmp)
    return layer_op(tensor)


def flatten(x):
    """Flatten tensor x"""
    return reshape(x, (-1,))


def stack(tensor_list, axis=0, **kwargs):
    """Stack a list of tensors
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/stack>`_
    """
    return Kops.stack(tensor_list, axis=axis)


def concatenate(tensor_list, axis=-1, target_shape=None, name=None):
    """
    Concatenates a list of numbers or tensor into a bigger tensor
    If the target shape is given, the output is reshaped to said shape
    """
    concatenated_tensor = Kops.concatenate(tensor_list, axis=axis)

    if target_shape is None:
        return concatenated_tensor
    return K.reshape(concatenated_tensor, target_shape)


def scatter_to_one(values, indices, output_shape):
    """
    Like scatter_nd initialized to one instead of zero
    see full `docs <https://www.tensorflow.org/api_docs/python/tf/scatter_nd>`_
    """
    ones = Kops.ones(output_shape)
    return Kops.scatter_update(ones, indices, values)


def swapaxes(tensor, source, destination):
    """
    Moves the axis of the tensor from source to destination, as in numpy.swapaxes.
    see full `docs <https://numpy.org/doc/stable/reference/generated/numpy.swapaxes.html>`_
    """
    rank = len(tensor.shape)
    indices = list(range(rank))
    if source < 0:
        source += rank
    if destination < 0:
        destination += rank

    indices[source], indices[destination] = indices[destination], indices[source]

    return Kops.transpose(tensor, indices)


def elu(x, alpha=1.0, **kwargs):
    new_layer = ELU(alpha=alpha, **kwargs)
    return new_layer(x)


def tensor_splitter(ishape, split_sizes, axis=2, name="splitter"):
    """
    Generates a Lambda layer to apply the split operation to a given tensor shape.
    This wrapper cannot split along the batch index (axis=0).

    Parameters
    ----------
        ishape: list(int)
            input shape of the tensor that will be split
        split_sizes: list(int)
            size of each chunk
        axis: int
            axis along which the split will be applied
        name: str
            name of the layer
    Returns
    -------
        sp_layer: layer
            a keras layer that applies the split operation upon call
    """
    if axis < 1:
        raise ValueError("tensor_splitter wrapper can only split along non-batch dimensions")

    # Check that we can indeed split this
    if ishape[axis] != np.sum(split_sizes):
        raise ValueError(
            f"Cannot split tensor of shape {ishape} along axis {axis} in chunks of {split_sizes}"
        )

    # Output shape of each split
    oshapes = []
    # Indices at which to put the splits
    # NB: tensorflow's split function would've taken the split_sizes directly
    # keras instead takes the index at where to split
    indices = []
    current_idx = 0

    for xsize in split_sizes:
        current_idx += xsize
        indices.append(current_idx)
        oshapes.append((*ishape[1:axis], xsize, *ishape[axis + 1 :]))

    sp_layer = keras_Lambda(
        lambda x: Kops.split(x, indices, axis=axis), output_shape=oshapes, name=name
    )
    return sp_layer
