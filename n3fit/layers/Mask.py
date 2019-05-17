import numpy as np
from backends import MetaLayer

class Mask(MetaLayer):
    """
    Mask is a layer used to apply both a boolean mask and a multiplier
    """
    def __init__(self, bool_mask, c = 1.0, **kwargs):
        self.output_dim = np.count_nonzero(bool_mask)
        self.mask = bool_mask
        self.c = c
        super(MetaLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        initializer = MetaLayer.init_constant(value = self.c)
        self.kernel = self.builder_helper("mask", (1,), initializer, trainable = False)
        super(Mask, self).build(input_shape)

    def compute_output_shape(self, input_shape):
        return (self.output_dim, None)

    def call(self, prediction_in):
        return self.boolean_mask(self.kernel*prediction_in, self.mask)


