import numpy as np
from n3fit.layers.Observable import Observable
from n3fit.backends import operations as op

import tensorflow as tf


class DIS(Observable):
    """
        The DIS class receives a list of active flavours and a fktable
        and prepares a layer that performs the convolution of said fktable with
        the incoming pdf.
    """

    def gen_basis(self, basis):
        """
            Receives a list of active flavours and generates a boolean mask

            # Arguments:
                - `basis`: list of active flavours
        """
        if basis is  None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def meta_call(self, pdf):
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
                    rank 1 tensor (ndata,)
        """
        # Do we need to add a splitting?
        if self.splitting is None:
            pdf_masked = op.boolean_mask(pdf, self.basis[0], axis = 2)
            pdfs_masked = [pdf_masked]*self.n_conv # TODO this can be written better
        else:
            pdfs_masked = []
            for partial_pdf, mask in zip(tf.split(pdf, self.splitting, axis=1), self.basis):
                pdfs_masked.append(op.boolean_mask(partial_pdf, mask, axis = 2))

        # Convolute with the fktable(s)
        results = []
        for fktable, pdf_masked in zip(self.fktables, pdfs_masked):
            res = op.tensor_product(pdf_masked, fktable, axes = [(1,2), (1,2)])
            results.append(res)

        # Apply the operation (if any)
        ret = self.operation(results)
        return ret
