from numpy import count_nonzero

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class Mask(MetaLayer):
    """
    This layers applies a boolean mask to a rank-1 input tensor.
    The mask admit a multiplier for all outputs which will be internally
    saved as a weight so it can be updated during trainig.

    Typical usage is to apply training/validation split masks
    or applying a multiplier to a given layer


    Parameters
    ----------
        bool_mask: np.array
            numpy array with the boolean mask to be applied
        c: float
            constant multiplier for every output
        axis: int
            axis in which to apply the mask. Currently,
            only the last axis gives the correct output shape
    """

    def __init__(self, bool_mask=None, c=None, axis=None, **kwargs):
        if bool_mask is None:
            self.mask = None
            self.last_dim = -1
        else:
            self.mask = op.numpy_to_tensor(bool_mask, dtype=bool)
            if len(bool_mask.shape) == 1:
                self.last_dim = count_nonzero(bool_mask)
            else:
                self.last_dim = count_nonzero(bool_mask[0, ...])
        self.c = c
        self.axis = axis
        super().__init__(**kwargs)

    def build(self, input_shape):
        if self.c is not None:
            initializer = MetaLayer.init_constant(value=self.c)
            self.kernel = self.builder_helper("mask", (1,), initializer, trainable=False)
        super(Mask, self).build(input_shape)

    def call(self, ret):
        if self.mask is not None:
            flat_res = op.boolean_mask(ret, self.mask, axis=self.axis)
            output_shape = [-1 if d is None else d for d in ret.get_shape()]
            # Make sure reshape will succeed: set the last dimension to the unmasked data length and before-last to
            # the number of replicas
            output_shape[-1] = self.last_dim
            output_shape[-2] = self.mask.shape[-2]
            ret = op.reshape(flat_res, shape=output_shape)
        if self.c is not None:
            ret = ret * self.kernel
        return ret
