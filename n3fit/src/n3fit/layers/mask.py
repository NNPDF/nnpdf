import numpy as np

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class Mask(MetaLayer):
    """
    This layers applies a boolean mask to an input tensor.
    The mask admit a multiplier for all outputs which will be internally
    saved as a weight so it can be updated during trainig.

    Typical usage is to apply training/validation split masks
    or applying a multiplier to a given layer


    Parameters
    ----------
        bool_mask: np.array of shape (n_replicas, n_features)
            numpy array with the boolean mask to be applied
        c: float
            constant multiplier for every output
    """

    def __init__(self, bool_mask=None, c=None, **kwargs):
        self._raw_mask = bool_mask
        self._flattened_indices = None
        if bool_mask is None:
            self.mask = None
            self.last_dim = -1
        else:
            self.mask = op.numpy_to_tensor(bool_mask, dtype=bool)
            self.last_dim = np.count_nonzero(bool_mask[0, ...])
        self.c = c
        self.masked_output_shape = None
        super().__init__(**kwargs)

    def build(self, input_shape):
        if self.c is not None:
            initializer = MetaLayer.init_constant(value=self.c)
            self.kernel = self.builder_helper("mask", (1,), initializer, trainable=False)
        # Make sure reshape will succeed: set the last dimension to the unmasked data length and before-last to
        # the number of replicas
        if self.mask is not None:

            # Prepare the indices to mask
            indices = np.where(self._raw_mask)
            # The batch dimension can be ignored
            nreps = self.mask.shape[-2]
            self._flattened_indices = np.ravel_multi_index(indices, self._raw_mask.shape)
            self.masked_output_shape = [-1 if d is None else d for d in input_shape]
            self.masked_output_shape[-1] = self.last_dim
            self.masked_output_shape[-2] = nreps
        super().build(input_shape)

    def call(self, ret):
        """
        Apply the mask to the input tensor, and multiply by the constant if present.

        Parameters
        ----------
            ret: Tensor of shape (batch_size, n_replicas, n_features)

        Returns
        -------
            Tensor of shape (batch_size, n_replicas, n_features)
        """
        if self.mask is not None:
            # TODO: JtH flattened indices are out of bounds: op.flatten(ret) = len(training) while self._flattened_indices.shape = training + val
            # import pdb; pdb.set_trace()
            ret = op.gather(op.flatten(ret), self._flattened_indices)
            ret = op.reshape(ret, self.masked_output_shape)
        if self.c is not None:
            ret = ret * self.kernel
        return ret
