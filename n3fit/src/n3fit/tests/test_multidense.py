import numpy as np
import tensorflow as tf
from tensorflow.keras import Sequential
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.layers import Dense

from n3fit.backends.keras_backend.multi_dense import MultiDense, MultiDropout
from n3fit.model_gen import generate_nn


def test_multidense():
    replicas = 2
    multi_dense_model = Sequential(
        [
            MultiDense(
                units=8,
                replica_seeds=[42, 43],
                replica_input=False,
                kernel_initializer=GlorotUniform(seed=0),
            ),
            MultiDense(units=4, replica_seeds=[52, 53], kernel_initializer=GlorotUniform(seed=100)),
        ]
    )
    single_models = [
        Sequential(
            [
                Dense(units=8, kernel_initializer=GlorotUniform(seed=42 + r)),
                Dense(units=4, kernel_initializer=GlorotUniform(seed=52 + r + 100)),
            ]
        )
        for r in range(replicas)
    ]

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
    dense_layers = []
    for r in range(2):
        dense_layer = Dense(units=2, kernel_initializer=GlorotUniform(seed=42 + r))
        dense_layer.build(input_shape=input_shape)
        dense_layers.append(dense_layer)
    stacked_weights = tf.stack([dense_layer.weights[0] for dense_layer in dense_layers], axis=0)

    multi_dense_layer = MultiDense(
        units=2,
        replica_seeds=[0, 1],
        replica_input=False,
        kernel_initializer=GlorotUniform(seed=42),
    )
    multi_dense_layer.build(input_shape=input_shape)

    multi_dense_weights = multi_dense_layer.weights[0].numpy()
    stacked_weights = stacked_weights.numpy()

    np.testing.assert_allclose(multi_dense_weights, stacked_weights)


def test_dropout():
    replicas = 100
    x_size = 10
    features = 1
    input_shape = (1, replicas, x_size, features)
    test_input = tf.ones(shape=input_shape)

    dropout_layer = MultiDropout(rate=0.5, seed=44)

    test_output = dropout_layer(test_input, training=True)

    # Check that for every replica the same x values are dropped
    zero_indices_first_replica = np.where(test_output.numpy()[0, 0, :, 0] == 0)
    nonzero_indices_first_replica = np.where(test_output.numpy()[0, 0, :, 0] != 0)

    assert np.all(test_output.numpy()[:, :, zero_indices_first_replica, :] == 0)
    assert np.all(test_output.numpy()[:, :, nonzero_indices_first_replica, :] != 0)
