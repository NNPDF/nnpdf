from n3fit.backends import MetaLayer
import tensorflow as tf


class Extrapolation(MetaLayer):
    """
        Nonlinear transformation of the x-grid, based on wether or not x is inside the data fk xgrids domain.
    """

    def __init__(
        self, smallxlim_scaled, **kwargs
    ):
        self.smallxlim_scaled = smallxlim_scaled
        super().__init__(**kwargs)

    def call(self, x):
        less_mask = tf.math.less(x, self.smallxlim_scaled)
        x0 = tf.keras.backend.ones_like(x) * self.smallxlim_scaled
        ret = tf.where(less_mask, x0, x)
        return ret
