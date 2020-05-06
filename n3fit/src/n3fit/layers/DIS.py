"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.
"""

import numpy as np
from n3fit.layers.Observable import Observable
from n3fit.backends import operations as op



class DIS(Observable):
    """
        The DIS class receives a list of active flavours and a fktable
        and prepares a layer that performs the convolution of said fktable with
        the incoming pdf.

        The fktable is expected to be rank 3 (ndata, xgrid, flavours)
        while the input pdf is also rank 3 where the first dimension is the batch dimension
        (1, xgrid, flavours)
    """

    def gen_mask(self, basis):
        """
            Receives a list of active flavours and generates a boolean mask tensor

            Parameters
            ----------
                basis: list(int)
                    list of active flavours
        """
        if basis is  None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def meta_call(self, pdf_raw):
        """
            Thiss function perform the fktable \otimes pdf convolution.

            Firs pass the input PDF through a mask to remove the unactive flavour, then transpose the PDF
            to have everything in the correct order and finally perform a tensorproduct contracting both pdf indices.

            Parameters
            ----------
                pdf:  backend tensor
                    rank 3 tensor (batch_size, xgrid, flavours)

            Returns
            -------
                result: backend tensor
                    rank 1 tensor (batchsize, ndata)
        """
        # DIS never needs splitting
        if self.splitting is not None:
            raise ValueError("DIS layer call with a dataset that needs more than one xgrids?")

        pdf = op.unbatch(pdf_raw)

        results = []

        # Separate the two possible paths this layer can take
        if self.many_masks:
            for mask, fktable in zip(self.all_masks, self.fktables):
                pdf_masked = op.boolean_mask(pdf, mask, axis=1)
                res = op.tensor_product(pdf_masked, fktable, axes = [(0,1), (2,1)])
                results.append(res)
        else:
            pdf_masked = op.boolean_mask(pdf, self.all_masks[0], axis=1)
            for mask, fktable in zip(self.all_masks, self.fktables):
                res = op.tensor_product(pdf_masked, fktable, axes = [(0,1), (2,1)])
                results.append(res)

        ret = self.operation(results)
        return op.batchit(ret)
