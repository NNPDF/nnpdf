from typing import List

import tensorflow as tf
from tensorflow.keras.initializers import Initializer
from tensorflow.keras.layers import Dense


class MultiDense(Dense):
    """
    Dense layer for multiple replicas at the same time.

    Inputs to this layer may contain multiple replicas, for the later layers.
    In this case, the `replica_input` argument should be set to `True`, which is the default.
    The input shape in this case is (batch_size, replicas, gridsize, features).
    For the first layer, there are no replicas yet, and so the `replica_input` argument
    should be set to `False`.
    The input shape in this case is (batch_size, gridsize, features).

    Weights are initialized using a `replica_seeds` list of seeds, and are identical to the
    weights of a list of single dense layers with the same `replica_seeds`.


    Example
    -------

    >>> from tensorflow.keras import Sequential
    >>> from tensorflow.keras.layers import Dense
    >>> from tensorflow.keras.initializers import GlorotUniform
    >>> import tensorflow as tf
    >>> replicas = 2
    >>> multi_dense_model = Sequential([
    >>>     MultiDense(units=8, replica_seeds=[42, 43], replica_input=False, kernel_initializer=GlorotUniform(seed=0)),
    >>>     MultiDense(units=4, replica_seeds=[52, 53], kernel_initializer=GlorotUniform(seed=0)),
    >>>     ])
    >>> single_models = [
    >>>     Sequential([
    >>>         Dense(units=8, kernel_initializer=GlorotUniform(seed=42 + r)),
    >>>         Dense(units=4, kernel_initializer=GlorotUniform(seed=52 + r)),
    >>>         ])
    >>>     for r in range(replicas)
    >>>     ]
    >>> gridsize, features = 100, 2
    >>> multi_dense_model.build(input_shape=(None, gridsize, features))
    >>> for single_model in single_models:
    >>>     single_model.build(input_shape=(None, gridsize, features))
    >>> test_input = tf.random.uniform(shape=(1, gridsize, features))
    >>> multi_dense_output = multi_dense_model(test_input)
    >>> single_dense_output = tf.stack([single_model(test_input) for single_model in single_models], axis=1)
    >>> tf.reduce_all(tf.equal(multi_dense_output, single_dense_output))

    Parameters
    ----------
    replica_seeds: List[int]
        List of seeds per replica for the kernel initializer.
    kernel_initializer: Initializer
        Initializer class for the kernel.
    replica_input: bool (default: True)
        Whether the input already contains multiple replicas.
    """

    def __init__(
        self,
        replica_seeds: List[int],
        kernel_initializer: Initializer,
        replica_input: bool = True,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.replicas = len(replica_seeds)
        self.replica_seeds = replica_seeds
        self.kernel_initializer = MultiInitializer(
            single_initializer=kernel_initializer, replica_seeds=replica_seeds
        )
        self.bias_initializer = MultiInitializer(
            single_initializer=self.bias_initializer, replica_seeds=replica_seeds
        )
        self.replica_input = replica_input

        self.matmul = None

    def build(self, input_shape):
        input_dim = input_shape[-1]
        self.kernel = self.add_weight(
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

        if self.replicas == 1:
            if self.replica_input:
                self.matmul = lambda inputs: tf.tensordot(inputs, self.kernel[0], [[-1], [0]])
            else:
                # Manually add replica dimension
                self.matmul = lambda inputs: tf.expand_dims(
                    tf.tensordot(inputs, self.kernel[0], [[-1], [0]]), axis=1
                )
        else:
            einrule = "brnf,rfg->brng" if self.replica_input else "bnf,rfg->brng"
            self.matmul = lambda inputs: tf.einsum(einrule, inputs, self.kernel)

    def call(self, inputs):
        """
        Compute output of shape (batch_size, replicas, gridsize, units).

        For the first layer, (self.replica_input is False), this is equivalent to
        applying each replica separately and concatenating along the last axis.
        If the input already contains multiple replica outputs, it is equivalent
        to applying each replica to its corresponding input.
        """
        if inputs.dtype.base_dtype != self._compute_dtype_object.base_dtype:
            inputs = tf.cast(inputs, dtype=self._compute_dtype_object)

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
        # Remove the replica axis from the input shape.
        if self.replica_input:
            input_shape = input_shape[:1] + input_shape[2:]

        output_shape = super().compute_output_shape(input_shape)

        # Add back the replica axis to the output shape.
        output_shape = output_shape[:1] + [self.replicas] + output_shape[1:]

        return output_shape

    def get_config(self):
        config = super().get_config()
        config.update({"replica_input": self.replica_input, "replica_seeds": self.replica_seeds})
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
    """

    def __init__(self, single_initializer: Initializer, replica_seeds: List[int]):
        self.initializer_class = type(single_initializer)
        self.initializer_config = single_initializer.get_config()
        self.base_seed = single_initializer.seed if hasattr(single_initializer, "seed") else None
        self.replica_seeds = replica_seeds

    def __call__(self, shape, dtype=None, **kwargs):
        shape = shape[1:]  # Remove the replica axis from the shape.
        per_replica_weights = []
        for replica_seed in self.replica_seeds:
            if self.base_seed is not None:
                self.initializer_config["seed"] = self.base_seed + replica_seed
            single_initializer = self.initializer_class.from_config(self.initializer_config)

            per_replica_weights.append(single_initializer(shape, dtype, **kwargs))

        return tf.stack(per_replica_weights, axis=0)
