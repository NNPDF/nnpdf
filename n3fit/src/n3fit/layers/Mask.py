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
        c: float
            constant multiplier for every output
        unbatch: bool
            whether to remove the first dimension (size=1) at the beginning
            if given sets batch_it to False
        axis: int
            axis in which to apply the mask
    """

    def __init__(self, bool_mask = None, c=None, unbatch=False, axis=None, **kwargs):
        self.output_dim = np.count_nonzero(bool_mask)
        if bool_mask is None:
            self.mask = None
        else:
            self.mask = op.numpy_to_tensor(bool_mask, dtype=bool)
        self.c = c
        self.axis = axis
        self.unbatch = unbatch
        super().__init__(**kwargs)

    def build(self, input_shape):
        if self.c is not None:
            initializer = MetaLayer.init_constant(value=self.c)
            self.kernel = self.builder_helper(
                "mask", (1,), initializer, trainable=False
            )
        super(Mask, self).build(input_shape)

    def meta_call(self, ret):
        if self.mask is not None:
            ret = op.boolean_mask(ret, self.mask, axis=self.axis)
        if self.c:
            ret = ret * self.kernel
        if self.unbatch:
            ret = op.unbatch(ret)
        return ret
