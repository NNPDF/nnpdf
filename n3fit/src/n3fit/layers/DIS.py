"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.
"""

import numpy as np

from n3fit.backends import operations as op

from .observable import Observable


class DIS(Observable):
    """
    The DIS class receives a list of active flavours and a fktable
    and prepares a layer that performs the convolution of said fktable with
    the incoming pdf.

    The fktable is expected to be rank 3 (ndata, xgrid, flavours)
    while the input pdf is rank 4 of shape (batch_size, replicas, xgrid, flavours)
    """

    def gen_mask(self, basis):
        """
        Receives a list of active flavours and generates a boolean mask tensor

        Parameters
        ----------
            basis: list(int)
                list of active flavours

        Returns
        -------
            mask: tensor
                rank 1 tensor (flavours)
        """
        if basis is None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def mask_fk(self, fk, mask):
        return op.einsum('fF, nFx -> nxf', mask, fk)

    def build(self, input_shape):
        super().build(input_shape)
        if self.num_replicas > 1:

            def compute_observable(pdf, masked_fk):
                return op.einsum('brxf, nxf -> brn', pdf, masked_fk)

        else:

            def compute_observable(pdf, masked_fk):
                return op.tensor_product(pdf, masked_fk, axes=[(2, 3), (1, 2)])  # brxf, nxf -> brn

        self.compute_observable = compute_observable
