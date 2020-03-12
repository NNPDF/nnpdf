import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations


class Mask(MetaLayer):
    """
    This layers applies a boolean mask to a rank-1 input tensor.
    By default returns a rank-2 tensor where the first dimension
    is a batch dimension of size 1

    Its most common use is the training/validation split.
    The mask admit a multiplier for all outputs

    Parameters
    ----------
        bool_mask: np.array
            numpy array with the boolean mask to be applied
        c: bool
            constant multiplier for every output
    """

    def __init__(self, bool_mask, c=1.0, batch_it = True, **kwargs):
        self.output_dim = np.count_nonzero(bool_mask)
        self.mask = bool_mask
        self.c = c
        self.batch_it = batch_it
        super(MetaLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        initializer = MetaLayer.init_constant(value=self.c)
        self.kernel = self.builder_helper("mask", (1,), initializer, trainable=False)
        super(Mask, self).build(input_shape)

    def call(self, prediction_in):
        ret = self.boolean_mask(self.kernel * prediction_in, self.mask)
        if self.batch_it:
            return operations.batchit(ret)
