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
    and the last dimension the number of replicas being fitted (1, xgrid, flavours, replicas)
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

    def call(self, pdf, A_to_idx):
        """
        This function perform the fktable \otimes pdf convolution.

        First pass the input PDF through a mask to remove the unactive flavors,
        then a tensor_product between the PDF and each fktable is performed
        finally the defined operation is applied to all the results

        Parameters
        ----------
            pdf:  list
                list of rank 4 tensors (batch_size, xgrid, flavours, replicas)
                for different A values
            A_to_idx: dict
                dictionary that maps A values to the index in the pdf

        Returns
        -------
            result: backend tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        # Select the correct A values from the pdf
        pdfs = super().call(pdf, A_to_idx)

        # DIS never needs splitting
        if self.splitting is not None:
            raise ValueError("DIS layer call with a dataset that needs more than one xgrid?")

        results = []
        # Separate the two possible paths this layer can take
        if self.many_masks:
            # TODO: not sure about this case?
            pdf = pdfs[0]
            for mask, fktable in zip(self.all_masks, self.fktables):
                pdf_masked = op.boolean_mask(pdf, mask, axis=2)
                res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                results.append(res)
        else:
            for pdf, fktable in zip(pdfs, self.fktables):
                pdf_masked = op.boolean_mask(pdf, self.all_masks[0], axis=2)
                res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                results.append(res)

        return self.operation(results)
