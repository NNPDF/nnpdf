"""
    Implementations of weight constraints for initializers
"""

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.constraints import MinMaxNorm


class MinMaxWeight(MinMaxNorm):
    """
    Small override to the MinMaxNorm Keras class to not look at the absolute value
    This version looks at the sum instead of at the norm
    """

    def __init__(self, min_value, max_value, **kwargs):
        super(MinMaxWeight, self).__init__(min_value=min_value, max_value=max_value, **kwargs)

    @tf.function
    def __call__(self, w):
        norms = K.sum(w, axis=self.axis, keepdims=True)
        desired = (
            self.rate * K.clip(norms, self.min_value, self.max_value) + (1 - self.rate) * norms
        )
        return w * desired / (K.epsilon() + norms)
