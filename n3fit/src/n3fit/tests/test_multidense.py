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
                replicas=replicas,
                seed=42,
                replica_input=False,
                initializer_class=GlorotUniform,
            ),
            MultiDense(units=4, replicas=replicas, seed=52, initializer_class=GlorotUniform),
        ]
    )
    single_models = [
        Sequential(
            [
                Dense(units=8, kernel_initializer=GlorotUniform(seed=42 + r)),
                Dense(units=4, kernel_initializer=GlorotUniform(seed=52 + r)),
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

    np.allclose(multi_dense_output, single_dense_output)
