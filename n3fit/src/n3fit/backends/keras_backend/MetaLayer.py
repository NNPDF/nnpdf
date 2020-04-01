"""
    The class MetaLayer is an extension of the backend Layer class
    with a number of methods and helpers to facilitate writing new custom layers
    in such a way that the new custom layer don't need to rely in anything backend-dependent

    In other words, if you want to implement a new layer and need functions not included here
    it is better to add a new method which is just a call to the relevant backend-dependent function
    For instance: np_to_tensor is just a call to K.constant
"""

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Layer
from tensorflow.keras.initializers import (
    Constant,
    RandomUniform,
    glorot_normal,
    glorot_uniform,
)

# Define in this dictionary new initializers as well as the arguments they accept (with default values if needed be)
initializers = {
    "random_uniform": (RandomUniform, {"minval": -0.5, "maxval": 0.5}),
    "glorot_uniform": (glorot_uniform, {}),
    "glorot_normal": (glorot_normal, {}),
}


def reinitialize_weight(weight, initializer, seed_skip=10):
    new_config = initializer.get_config()
    if "seed" in new_config:
        new_seed = initializer.seed + 10
        new_config["seed"] += new_seed
        initializer.seed = new_seed
    new_init = initializer.from_config(new_config)
    new_weights = new_init(weight.shape)
    weight.assign(new_weights)


class MetaLayer(Layer):
    """
    This metalayer function must contain all backend-dependent functions

    In order to write a custom Keras layer you usually need to override:
        - __init__
        - build
        - call
        - output_shape
    """

    initializers = initializers
    weight_inits = []

    # Building function
    def builder_helper(
        self, name, kernel_shape, initializer, trainable=True, constraint=None
    ):
        """
        Creates a kernel that should be saved as an attribute of the caller class
        name: name of the kernel
        shape: tuple with its shape
        initializer: one of the initializers from this class (actually, any keras initializer)
        trainable: if it is
        constraint: one of the constraints from this class (actually, any keras constraints)
        """
        kernel = self.add_weight(
            name=name,
            shape=kernel_shape,
            initializer=initializer,
            trainable=trainable,
            constraint=constraint,
        )
        if trainable:
            self.weight_inits.append((kernel, initializer))
        return kernel

    def reinitialize(self):
        """ Reinitialize all weights and kernels with an initializer """
        # First check the default keras-tf ones
        if hasattr(self, "kernel_initializer"):
            kern = self.kernel
            kern_init = self.kernel_initializer
            bias = self.bias
            bias_init = self.bias_initializer
            reinitialize_weight(kern, kern_init)
            reinitialize_weight(bias, bias_init)
        # Now check the ones generated with the builder helper
        for w, init in self.weight_inits:
            reinitialize_weight(w, init)

    # Implemented initializers
    @staticmethod
    def init_constant(value):
        return Constant(value=value)

    @staticmethod
    def select_initializer(ini_name, seed=None, **kwargs):
        """
        Selects one of the initializers (which does initialize, i.e., not constant)
        All of them should accept seed
        """
        try:
            ini_tuple = initializers[ini_name]
        except KeyError as e:
            raise NotImplementedError(
                f"[MetaLayer.select_initializer] initializer not implemented: {ini_name}"
            ) from e

        ini_class = ini_tuple[0]
        ini_args = ini_tuple[1]
        ini_args["seed"] = seed

        for key, value in kwargs.items():
            if key in ini_args.keys():
                ini_args[key] = value
        return ini_class(**ini_args)

    # Make numpy array into a tensor
    def np_to_tensor(self, np_array, **kwargs):
        """
        Given a numpy array, returns a constant tensor
        """
        return K.constant(np_array, **kwargs)

    # Common tensor operations
    def tensor_ones(self, shape, **kwargs):
        """
        Generates a tensor of ones of the given shape
        """
        return K.ones(shape, **kwargs)

    def tensor_ones_like(self, tensor, **kwargs):
        """
        Generates a tensor of ones of the same shape as the input tensor
        """
        return K.ones_like(tensor, **kwargs)

    def many_replication(self, grid, replications, axis=0, **kwargs):
        """
            Generates a tensor with one extra dimension:
                a repetition of "grid" n times along the given axis
            from keras documentation:
            If x has shape (s1, s2, s3) and axis is 1, the output will have shape (s1, s2 * rep, s3)
        """
        return K.repeat_elements(grid, rep=replications, axis=axis, **kwargs)

    def sum(self, tensor, axis=None, **kwargs):
        """
        Computes the sum of the elements of the tensor
        """
        return K.sum(tensor, axis=axis, **kwargs)

    def tensor_product(self, tensor_x, tensor_y, axes, **kwargs):
        """
        Computes the tensordot product between tensor_x and tensor_y
        """
        return tf.tensordot(tensor_x, tensor_y, axes=axes, **kwargs)

    def transpose(self, tensor, **kwargs):
        """
        Transpose a tensor
        """
        return K.transpose(tensor, **kwargs)

    def boolean_mask(self, tensor, mask, axis=None, **kwargs):
        """
        Applies boolean mask to a tensor
        """
        return tf.boolean_mask(tensor, mask, axis=axis, **kwargs)

    def concatenate(self, tensor_list, axis=-1, target_shape=None):
        """
        Concatenates a list of numbers or tenosr into a bigger tensor
        If the target shape is given, the output is reshaped to said shape
        """
        concatenated_tensor = K.concatenate(tensor_list, axis=axis)
        if target_shape:
            return K.reshape(concatenated_tensor, target_shape)
        else:
            return concatenated_tensor

    def flatten(self, x):
        """ Flatten tensor x """
        return tf.reshape(x, (-1,))

    def permute_dimensions(self, tensor, permutation, **kwargs):
        """
        Receives a tensor and a tuple and permutes the axes of the tensor according to it.
        i.e.
        if permutation = (1,0,2)
        does the permutation: axis_0 -> axis_1, axis_1 -> axis_0, axis_2 -> axis_2
        """
        return K.permute_dimensions(tensor, permutation, **kwargs)
