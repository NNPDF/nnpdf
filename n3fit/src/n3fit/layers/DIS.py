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
    while the input pdf is rank 4 where the first dimension is the batch dimension
    and the last dimension the number of replicas being fitted (1, replicas, xgrid, flavours)
    """

    def gen_mask(self, basis):
        """
        Receives a list of active flavours and generates a boolean mask tensor

        Parameters
        ----------
            basis: list(int)
                list of active flavours
        """
        if basis is None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def call(self, pdf):
        """
        This function perform the fktable \otimes pdf convolution.

        First pass the input PDF through a mask to remove the unactive flavors,
        then a tensor_product between the PDF and each fktable is performed
        finally the defined operation is applied to all the results

        Parameters
        ----------
            pdf:  backend tensor
                rank 4 tensor (batch_size, replicas, xgrid, flavours)

        Returns
        -------
            result: backend tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        # DIS never needs splitting
        if self.splitting is not None:
            raise ValueError("DIS layer call with a dataset that needs more than one xgrid?")

        results = []
        # Separate the two possible paths this layer can take
        if self.many_masks:
            for mask, fktable in zip(self.all_masks, self.fktables):
                pdf_masked = op.boolean_mask(pdf, mask, axis=3)
                res = op.tensor_product(pdf_masked, fktable, axes=[(2, 3), (2, 1)])
                results.append(res)
        else:
            pdf_masked = op.boolean_mask(pdf, self.all_masks[0], axis=3)
            for fktable in self.fktables:
                res = op.tensor_product(pdf_masked, fktable, axes=[(2, 3), (2, 1)])
                results.append(res)

        return self.operation(results)
