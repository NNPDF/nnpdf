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
<<<<<<< HEAD
        c: bool
            constant multiplier for every output
    """

    def __init__(self, bool_mask, c=1.0, batch_it = True, **kwargs):
=======
        c: float
            constant multiplier for every output
        unbatch: bool
            whether to remove the first dimension (size=1) at the beginning
            if given sets batch_it to False
        axis: int
            axis in which to apply the mask
    """

    def __init__(self, bool_mask, c=None, unbatch=False, axis=None, **kwargs):
>>>>>>> master
        self.output_dim = np.count_nonzero(bool_mask)
        self.mask = op.numpy_to_tensor(bool_mask)
        self.c = c
<<<<<<< HEAD
        self.batch_it = batch_it
=======
        self.axis = axis
        self.unbatch = unbatch
>>>>>>> master
        super().__init__(**kwargs)

    def build(self, input_shape):
        if self.c:
            initializer = MetaLayer.init_constant(value=self.c)
            self.kernel = self.builder_helper(
                "mask", (1,), initializer, trainable=False
            )
        super(Mask, self).build(input_shape)

<<<<<<< HEAD
    def call(self, prediction_in):
        ret = op.boolean_mask(self.kernel * prediction_in, self.mask)
        if self.batch_it:
            return op.batchit(ret)
        else:
            return ret
=======
    def meta_call(self, prediction_in):
        ret = op.boolean_mask(prediction_in, self.mask, axis=self.axis)
        if self.c:
            ret = ret * self.kernel
        if self.unbatch:
            ret = op.unbatch(ret)
        return ret
>>>>>>> master
