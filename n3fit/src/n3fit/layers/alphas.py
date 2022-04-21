import tensorflow as tf
from tensorflow.keras.layers import Layer
import tensorflow_probability as tfp

min_alphas = 0.114
max_alphas = 0.122

random_initializer=tf.keras.initializers.RandomUniform(
    minval=min_alphas, 
    maxval=max_alphas,
)

class AlphasLayer(Layer):
    def __init__(self, seed, **kwargs):
        random_initializer=tf.keras.initializers.RandomUniform(
            minval=min_alphas, 
            maxval=max_alphas,
            seed=seed,
        )
        self.alphas = tf.Variable(
            initial_value=random_initializer(shape=(), dtype=tf.float32),
            trainable=True,
            name="alphas",
            dtype=tf.float32,
            constraint=lambda z: tf.clip_by_value(z, min_alphas, max_alphas)
        )
        super().__init__(**kwargs)

    def call(self, list_alphas_results, **kwargs):
        out = tfp.math.interp_regular_1d_grid(
            self.alphas,
            min_alphas,
            max_alphas,
            list_alphas_results,
            axis=0
        )
        return out
