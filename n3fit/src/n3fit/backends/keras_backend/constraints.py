"""
    Implementations of weight constraints for initializers
"""

from keras import backend as K
from keras.constraints import MinMaxNorm

from . import operations as ops


class MinMaxWeight(MinMaxNorm):
    """
    Small override to the MinMaxNorm Keras class to not look at the absolute value
    This version looks at the sum instead of at the norm
    """

    def __init__(self, min_value, max_value, **kwargs):
        super().__init__(min_value=min_value, max_value=max_value, axis=1, **kwargs)

    def __call__(self, w):
        norms = ops.sum(w, axis=self.axis, keepdims=True)
        desired = (
            self.rate * ops.clip(norms, self.min_value, self.max_value) + (1 - self.rate) * norms
        )
        return w * desired / (K.epsilon() + norms)
