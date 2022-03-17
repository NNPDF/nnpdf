import tensorflow as tf
# ALPHAS = tf.Variable(
#     initial_value=0.118,
#     trainable=True,
#     name="alphas",
#     dtype=tf.float32,
# )

from tensorflow.keras.layers import Layer
class AlphasLayer(Layer):
    def __init__(self, **kwargs):
        self.alphas = tf.Variable(
        initial_value=0.118,
        trainable=True,
        name="alphas",
        dtype=tf.float32,
    )
        super().__init__(**kwargs)

    def call(self, list_alphas_results, **kwargs):
        import tensorflow_probability as tfp
        out = tfp.math.interp_regular_1d_grid(
            self.alphas,
            0.116,
            0.120,
            list_alphas_results,
            fill_value="extrapolate",
            axis=0
        )
        return out