"""
The class MetaLayer is an extension of the backend Layer class
with a number of methods and helpers to facilitate writing new custom layers
in such a way that the new custom layer don't need to rely in anything backend-dependent

In other words, if you want to implement a new layer and need functions not included here
it is better to add a new method which is just a call to the relevant backend-dependent function
For instance: np_to_tensor is just a call to K.constant
"""

import math

from keras import random
from keras.initializers import Constant, RandomUniform, VarianceScaling, glorot_uniform
from keras.layers import Layer


class GammaVarianceScaling(VarianceScaling):
    """``VarianceScaling`` with a tunable exponent ``gamma`` on the variance.

    keras' ``VarianceScaling`` draws weights with variance ``scale / fan`` (standard
    deviation ``sqrt(scale / fan)``). This variant raises that variance to the power
    ``gamma``:

        variance = (scale / fan) ** gamma,   std = (scale / fan) ** (gamma / 2).

    ``gamma = 1`` reproduces ``VarianceScaling`` (hence ``glorot_normal``) exactly;
    ``gamma > 1`` makes the initialisation narrower, ``gamma < 1`` wider. ``fan`` is
    ``fan_in``, ``fan_out`` or their average, per ``mode`` (``fan_avg`` for
    ``glorot_normal``).
    """

    # keras' correction so a truncated normal has the requested std after truncation.
    _TRUNCATED_CORRECTION = 0.87962566103423978

    def __init__(
        self, gamma=1.0, scale=1.0, mode="fan_in", distribution="truncated_normal", seed=None
    ):
        super().__init__(scale=scale, mode=mode, distribution=distribution, seed=seed)
        self.gamma = gamma

    @staticmethod
    def _compute_fans(shape):
        """fan_in, fan_out for a weight of the given shape (matches keras)."""
        shape = tuple(shape)
        if len(shape) < 1:
            fan_in = fan_out = 1
        elif len(shape) == 1:
            fan_in = fan_out = shape[0]
        elif len(shape) == 2:
            fan_in, fan_out = shape
        else:
            receptive_field_size = 1
            for dim in shape[:-2]:
                receptive_field_size *= dim
            fan_in = shape[-2] * receptive_field_size
            fan_out = shape[-1] * receptive_field_size
        return float(fan_in), float(fan_out)

    def __call__(self, shape, dtype=None):
        scale = self.scale
        fan_in, fan_out = self._compute_fans(shape)
        if self.mode == "fan_in":
            scale /= max(1.0, fan_in)
        elif self.mode == "fan_out":
            scale /= max(1.0, fan_out)
        else:
            scale /= max(1.0, (fan_in + fan_out) / 2.0)
        # `scale` is now the post-division variance scale/fan; keras would take
        # std = sqrt(scale). Raise the *variance* to gamma -> std = scale**(gamma/2).
        # gamma=1 gives sqrt(scale) (standard glorot).
        std = scale ** (self.gamma / 2.0)
        if self.distribution == "truncated_normal":
            return random.truncated_normal(
                shape,
                mean=0.0,
                stddev=std / self._TRUNCATED_CORRECTION,
                dtype=dtype,
                seed=self.seed,
            )
        elif self.distribution == "untruncated_normal":
            return random.normal(shape, mean=0.0, stddev=std, dtype=dtype, seed=self.seed)
        else:  # uniform: keras uses limit = sqrt(3 * variance) = sqrt(3) * std
            limit = math.sqrt(3.0) * std
            return random.uniform(shape, minval=-limit, maxval=limit, dtype=dtype, seed=self.seed)

    def get_config(self):
        return {**super().get_config(), "gamma": self.gamma}


# Define in this dictionary new initializers as well as the arguments they accept (with default values if needed be)
initializers = {
    "random_uniform": (RandomUniform, {"minval": -0.5, "maxval": 0.5}),
    "glorot_uniform": (glorot_uniform, {}),
    # glorot_normal expressed via GammaVarianceScaling so its width is tunable through
    # `scale` (variance multiplier) and `gamma` (exponent on the variance:
    # variance = (scale/fan)**gamma). scale=1.0, gamma=1.0 reproduces keras'
    # glorot_normal exactly; gamma>1 narrower, gamma<1 wider.
    "glorot_normal": (
        GammaVarianceScaling,
        {"scale": 1.0, "gamma": 1.0, "mode": "fan_avg", "distribution": "untruncated_normal"},
    ),
}


class MetaLayer(Layer):
    """
    This metalayer function must contain all backend-dependent functions

    In order to write a custom Keras layer you usually need to override:
        - __init__
        - call
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
        main_name = f"{self.name}/{weight_name}"
        for weight in self.weights:
            if weight.name in (f"{main_name}:{internal_count}", main_name, weight_name):
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
        # Copy so per-call overrides (seed, scale, ...) don't leak into the shared defaults
        ini_args = dict(ini_tuple[1])
        ini_args["seed"] = seed

        for key, value in kwargs.items():
            if key in ini_args:
                ini_args[key] = value
        return ini_class(**ini_args)
