import tensorflow as tf
from tensorflow.keras.layers import Layer

class AlphasLayer(Layer):
    def __init__(self, **kwargs):
        self.alphas = tf.Variable(
        initial_value=0.118,
        trainable=True,
        name="alphas",
        dtype=tf.float32,
        constraint=lambda z: tf.clip_by_value(z, 0.114, 0.122)
    )
        super().__init__(**kwargs)

    def call(self, list_alphas_results, **kwargs):
        import tensorflow_probability as tfp
        out = tfp.math.interp_regular_1d_grid(
            self.alphas,
            0.114,
            0.122,
            list_alphas_results,
            axis=0
        )
        return out