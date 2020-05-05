import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


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
        batch_it: bool
            whether to add a size=1 dimension at the beginning
        unbatch: bool
            whether to remove the first dimension (size=1) at the beginning
            if given sets batch_it to False
        axis: int
            axis in which to apply the mask
    """

    def __init__(self, bool_mask, c=1.0, batch_it = True, unbatch = False, axis = None, **kwargs):
        self.output_dim = np.count_nonzero(bool_mask)
        self.mask = bool_mask
        self.c = c
        self.batch_it = batch_it and not unbatch
        self.axis = axis
        self.unbatch = unbatch
        super().__init__(**kwargs)

    def build(self, input_shape):
        initializer = MetaLayer.init_constant(value=self.c)
        self.kernel = self.builder_helper("mask", (1,), initializer, trainable=False)
        super(Mask, self).build(input_shape)

    def meta_call(self, prediction_in):
        ret = op.boolean_mask(self.kernel * prediction_in, self.mask, axis = self.axis)
        if self.batch_it:
            ret = op.batchit(ret)
        if self.unbatch:
            ret = op.unbatch(ret)
        return ret
