import numpy as np
import tensorflow as tf
from tensorflow.keras import Sequential
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.layers import Dense

from n3fit.backends.keras_backend.multi_dense import MultiDense


def test_multidense():
    replicas = 2
    multi_dense_model = Sequential(
        [
            MultiDense(
                units=8,
                replica_seeds=[42, 43],
                is_first_layer=True,
                kernel_initializer=GlorotUniform(seed=5),
            ),
            MultiDense(units=4, replica_seeds=[52, 53], kernel_initializer=GlorotUniform(seed=100)),
        ]
    )
    single_models = []
    for r in range(replicas):
        single_models.append(
            Sequential(
                [
                    Dense(units=8, kernel_initializer=GlorotUniform(seed=42 + r + 5)),
                    Dense(units=4, kernel_initializer=GlorotUniform(seed=52 + r + 100)),
                ]
            )
        )

    gridsize, features = 100, 3
    multi_dense_model.build(input_shape=(None, gridsize, features))
    for single_model in single_models:
        single_model.build(input_shape=(None, gridsize, features))

    test_input = tf.random.uniform(shape=(1, gridsize, features))
    multi_dense_output = multi_dense_model(test_input)
    single_dense_output = tf.stack(
        [single_model(test_input) for single_model in single_models], axis=1
    )

    np.testing.assert_allclose(multi_dense_output, single_dense_output, atol=1e-6, rtol=1e-4)


def test_initializers():
    input_shape = (None, 3, 1)
    dense_weights = []
    for r in range(2):
        dense_layer = Dense(units=2, kernel_initializer=GlorotUniform(seed=42 + r))
        dense_layer.build(input_shape=input_shape)
        try:
            dense_weights.append(dense_layer.weights[0].value.numpy())
        except AttributeError:
            # In tensorflow < 2.16, value was a function
            dense_weights.append(dense_layer.weights[0].value().numpy())

    stacked_weights = np.stack(dense_weights, axis=0)

    multi_dense_layer = MultiDense(
        units=2,
        replica_seeds=[0, 1],
        is_first_layer=True,
        kernel_initializer=GlorotUniform(seed=42),
    )
    multi_dense_layer.build(input_shape=input_shape)

    multi_dense_weights = multi_dense_layer.weights[0].numpy()

    np.testing.assert_allclose(multi_dense_weights, stacked_weights)
