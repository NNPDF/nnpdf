"""
    The class MetaLayer is an extension of the backend Layer class
    with a number of methods and helpers to facilitate writing new custom layers
    in such a way that the new custom layer don't need to rely in anything backend-dependent

    In other words, if you want to implement a new layer and need functions not included here
    it is better to add a new method which is just a call to the relevant backend-dependent function
    For instance: np_to_tensor is just a call to K.constant
"""

from tensorflow.keras.initializers import Constant, RandomUniform, glorot_normal, glorot_uniform
from tensorflow.keras.layers import Layer

# Define in this dictionary new initializers as well as the arguments they accept (with default values if needed be)
initializers = {
    "random_uniform": (RandomUniform, {"minval": -0.5, "maxval": 0.5}),
    "glorot_uniform": (glorot_uniform, {}),
    "glorot_normal": (glorot_normal, {}),
}


class MetaLayer(Layer):
    """
    This metalayer function must contain all backend-dependent functions

    In order to write a custom Keras layer you usually need to override:
        - __init__
        - meta_call
    """

    initializers = initializers
    weight_inits = []

    # Building function
    def builder_helper(self, name, kernel_shape, initializer, trainable=True, constraint=None):
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

    def get_weight_by_name(self, weight_name, internal_count=0):
        """
        Returns a weight of the layer by name, returns None if the layer does not include
        the named weight.

        Note that internally weights of a layer are prefaced by the name of the layer, this
        should not be added to the input of this function. i.e., if the internal name is
        "layer/weight:0", the argument to this method should be just "weight".

        Parameters
        ----------
            weight_name: str
                Name of the weight
        """
        check_name = f"{self.name}/{weight_name}:{internal_count}"
        for weight in self.weights:
            if weight.name == check_name:
                return weight
        return None

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
