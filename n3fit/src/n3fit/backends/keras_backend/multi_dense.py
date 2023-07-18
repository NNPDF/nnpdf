import tensorflow as tf
from tensorflow.keras.initializers import Initializer
from tensorflow.keras.layers import Dense


class MultiDense(Dense):
    """
    Dense layer for multiple replicas at the same time.

    The kernel initializer is set using the custom arguments `initializer_class` and
    `seed`. The `seed` is incremented by 1 for each replica.

    Inputs to this layer may contain multiple replicas, for the later layers.
    In this case, the `replica_input` argument should be set to `True`, which is the default.
    For the first layer, there are no replicas yet, and so the `replica_input` argument
    should be set to `False`.

    Example
    -------

    >>> from tensorflow.keras import Sequential
    >>> from tensorflow.keras.layers import Dense
    >>> from tensorflow.keras.initializers import GlorotUniform
    >>> import tensorflow as tf
    >>> replicas = 2
    >>> multi_dense_model = Sequential([
    >>>     MultiDense(units=8, replicas=replicas, seed=42, replica_input=False, initializer_class=GlorotUniform),
    >>>     MultiDense(units=4, replicas=replicas, seed=52, initializer_class=GlorotUniform),
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
    >>> single_dense_output = tf.stack([single_model(test_input) for single_model in single_models], axis=-1)
    >>> tf.reduce_all(tf.equal(multi_dense_output, single_dense_output))

    Parameters
    ----------
    replicas: int
        Number of replicas.
    seed: int
        Seed for the random number generator.
    initializer_class: Initializer
        Initializer class for the kernel.
    replica_input: bool (default: True)
        Whether the input already contains multiple replicas.
    """
    def __init__(
        self,
        replicas: int,
        seed: int,
        initializer_class: Initializer,
        replica_input: bool = True,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.replicas = replicas
        self.seed = seed
        self.initializer_class = initializer_class
        self.replica_input = replica_input

    def build(self, input_shape):
        """
        Build weight matrix of shape (input_dim, units, replicas).
        Weights are initialized on a per-replica basis, with incrementing seed.
        """
        if self.replica_input:
            input_shape = input_shape[:-1]

        replica_kernels = []
        replica_biases = []
        for r in range(self.replicas):
            self.kernel_initializer = self.initializer_class(self.seed + r)
            super().build(input_shape)
            replica_kernels.append(self.kernel)
            replica_biases.append(self.bias)

        self.kernel = tf.Variable(tf.stack(replica_kernels, axis=-1))
        if self.use_bias:
            self.bias = tf.Variable(tf.stack(replica_biases, axis=-1))

        if self.replica_input:
            self.input_spec.axes = {-2: self.input_spec.axes[-1], -1: self.replicas}

    def call(self, inputs):
        """
        Compute output of shape (batch_size, gridsize, units, replicas).

        For the first layer, (self.replica_input is False), this is equivalent to
        applying each replica separately and concatenating along the last axis.
        If the input already contains multiple replica outputs, it is equivalent
        to applying each replica to its corresponding input.
        """
        if inputs.dtype.base_dtype != self._compute_dtype_object.base_dtype:
            inputs = tf.cast(inputs, dtype=self._compute_dtype_object)

        input_axes = 'bnfr' if self.replica_input else 'bnf'
        einrule = input_axes + ',fgr->bngr'
        outputs = tf.einsum(einrule, inputs, self.kernel)

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
        if self.replica_input:
            input_shape = input_shape[:-1]
        output_shape = super().compute_output_shape(input_shape)
        output_shape = output_shape.concatenate(self.replicas)
        return output_shape

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "replicas": self.replicas,
                "replica_input": self.replica_input,
                "seed": self.seed,
            }
        )
        return config
