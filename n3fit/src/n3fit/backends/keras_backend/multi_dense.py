"""
    Extend the ``Dense`` layer from Keras to act on an arbitrary number of replicas.
    This extension provides a performance improvement with respect to the original
    Dense layer from Keras even in the single replica case.
"""  # Tested last: Feb 2024

import tensorflow as tf
from tensorflow.keras.initializers import Initializer
from tensorflow.keras.layers import Dense

# Note for developers:
# This class plays with fire as it exploits the internals of Keras
# In particular, when saving the weights of the model, the current version (3.2)
# will rely on the existence of Dense being formed by two weights: a kernel and a bias
# the kernel variable is saved in the `_kernel` attribute while the bias  in `bias`.
# These should, in addition, correspond to weights "0" and "1" (and are saved as such).
# but this is not a public interface and so it could change at any point.
# In addition, from 3.2 the only accepted filename is `<name>.weight.h5`


class MultiDense(Dense):
    """
    Dense layer for multiple replicas at the same time.

    For the first layer in the network, (for which ``is_first_layer`` should be set to True),
    the input shape is (batch_size, gridsize, features), still without a replica axis.
    In this case this layer acts as a stack of single dense layers,
    with their own kernel and bias, acting on the same input.

    For subsequent layers, the input already contains multiple replicas, and the shape is
    (batch_size, replicas, gridsize, features).
    In this case, the input for each replica is multiplied by its own slice of the kernel.

    Weights are initialized using a `replica_seeds` list of seeds, and are identical to the
    weights of a list of single dense layers with the same `replica_seeds`.

    Parameters
    ----------
    replica_seeds: List[int]
        List of seeds per replica for the kernel initializer.
    kernel_initializer: Initializer
        Initializer class for the kernel.
    is_first_layer: bool (default: False)
        Whether this is the first MultiDense layer in the network, and so the input shape
        does not contain a replica axis.
    base_seed: int (default: 0)
        Base seed for the single replica initializer to which the replica seeds are added.
    """

    def __init__(
        self,
        replica_seeds: list[int],
        kernel_initializer: Initializer,
        is_first_layer: bool = False,
        base_seed: int = 0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.replicas = len(replica_seeds)
        self.replica_seeds = replica_seeds
        self.kernel_initializer = MultiInitializer(
            single_initializer=kernel_initializer, replica_seeds=replica_seeds, base_seed=base_seed
        )
        self.bias_initializer = MultiInitializer(
            single_initializer=self.bias_initializer,
            replica_seeds=replica_seeds,
            base_seed=base_seed,
        )
        self.is_first_layer = is_first_layer

        # Definition of the convolution between the input of the layer and the kernel parameters
        # it is defined during the build stage.
        self.matmul = None

        # See note above
        self._kernel = None
        self.bias = None

    def build(self, input_shape):
        input_dim = input_shape[-1]
        self._kernel = self.add_weight(
            name="kernel",
            shape=(self.replicas, input_dim, self.units),
            initializer=self.kernel_initializer,
            regularizer=self.kernel_regularizer,
            constraint=self.kernel_constraint,
        )
        if self.use_bias:
            self.bias = self.add_weight(
                name="bias",
                shape=(self.replicas, 1, self.units),
                initializer=self.bias_initializer,
                regularizer=self.bias_regularizer,
                constraint=self.bias_constraint,
            )
        else:
            self.bias = None
        self.input_spec.axes = {-1: input_dim}
        self.built = True

        # Using tensordot here for numerical stability with 4.0 fits
        # TODO: benchmark against the replica-agnostic einsum below and make that default
        # see https://github.com/NNPDF/nnpdf/pull/1905#discussion_r1489344081
        if self.replicas == 1:
            matmul = lambda inputs: tf.tensordot(inputs, self._kernel[0], [[-1], [0]])
            if self.is_first_layer:
                # Manually add replica dimension
                self.matmul = lambda x: tf.expand_dims(matmul(x), axis=1)
            else:
                self.matmul = matmul
        else:
            einrule = "bnf,rfg->brng" if self.is_first_layer else "brnf,rfg->brng"
            self.matmul = lambda inputs: tf.einsum(einrule, inputs, self._kernel)

    def call(self, inputs):
        """
        Compute output of shape (batch_size, replicas, gridsize, units).

        For the first layer, this is equivalent to
        applying each replica separately and concatenating along the last axis.
        If the input already contains multiple replica outputs, it is equivalent
        to applying each replica to its corresponding input.
        """
        # cast always
        inputs = tf.cast(inputs, dtype=self.compute_dtype)

        outputs = self.matmul(inputs)

        # Reshape the output back to the original ndim of the input.
        if not tf.executing_eagerly():
            output_shape = self.compute_output_shape(inputs.shape.as_list())
            outputs.set_shape(output_shape)

        if self.use_bias:
            outputs = outputs + self.bias

        if self.activation is not None:
            outputs = self.activation(outputs)

        return outputs

    def compute_output_shape(self, input_shape):
        if not self.is_first_layer:
            # Remove the replica axis from the input shape.
            input_shape = input_shape[:1] + input_shape[2:]

        output_shape = super().compute_output_shape(input_shape)

        # Add back the replica axis to the output shape.
        output_shape = output_shape[:1] + (self.replicas,) + output_shape[1:]

        return output_shape

    def get_config(self):
        config = super().get_config()
        config.update({"is_first_layer": self.is_first_layer, "replica_seeds": self.replica_seeds})
        return config


class MultiInitializer(Initializer):
    """
    Multi replica initializer that exactly replicates a stack of single replica initializers.

    Weights are stacked on the first axis, and per replica seeds are added to a base seed of the
    given single replica initializer.

    Parameters
    ----------
        single_initializer: Initializer
            Initializer class for the kernel.
        replica_seeds: List[int]
            List of seeds per replica for the kernel initializer.
        base_seed: int
            Base seed for the single replica initializer to which the replica seeds are added.
    """

    def __init__(self, single_initializer: Initializer, replica_seeds: list[int], base_seed: int):
        self.initializer_class = type(single_initializer)
        self.initializer_config = single_initializer.get_config()
        self.base_seed = base_seed
        self.replica_seeds = replica_seeds

    def __call__(self, shape, dtype=None, **kwargs):
        shape = shape[1:]  # Remove the replica axis from the shape.
        per_replica_weights = []
        for replica_seed in self.replica_seeds:
            if "seed" in self.initializer_config:
                self.initializer_config["seed"] = int(self.base_seed + replica_seed)
            single_initializer = self.initializer_class.from_config(self.initializer_config)

            per_replica_weights.append(single_initializer(shape, dtype, **kwargs))

        return tf.stack(per_replica_weights, axis=0)
