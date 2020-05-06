import numpy as np
from n3fit.layers.Observable import Observable
from n3fit.backends import operations as op

import tensorflow as tf


class DY(Observable):
    """
    Computes the convolution of two PDFs (the same one twice) and one fktable
    """

    def gen_basis(self, basis):
        if basis is None:
            basis_mask = np.ones((self.nfl, self.nfl), dtype=bool)
        else:
            basis_mask = np.zeros((self.nfl, self.nfl), dtype=bool)
            for i, j in basis.reshape(-1, 2):
                basis_mask[i, j] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def meta_call(self, pdf_raw):
        """
        This function perform the fktable \otimes pdf \otimes pdf convolution.

        First uses the basis of active combinations to generate a luminosity tensor
        with only some flavours active.

        The concatenate function returns a rank-3 tensor (combination_index, xgrid, xgrid)
        which can in turn be contracted with the rank-4 fktable.

        Parameters
        ----------
            pdf_in: tensor
                rank 3 tensor (batchsize, xgrid, flavours)

        Returns
        -------
            results: tensor
                rank 1 tensor (ndata,)
        """
        # Do we need to add a splitting?
        if self.splitting is None:
            pdf_x_pdf = op.pdf_masked_convolution(pdf_raw, self.basis[0])
            luminosities = [pdf_x_pdf]*self.n_conv
        else:
            luminosities = []
            for pdf, mask in zip(tf.split(pdf_raw, self.splitting, axis=1), self.basis):
                luminosities.append(op.pdf_masked_convolution(pdf, mask))

        results = []
        for fktable, luminosity in zip(self.fktables, luminosities):
            res = op.tensor_product(fktable, luminosity, axes=3)
            results.append(res)

        # Apply the operation (if any)
        ret = self.operation(results)
        return op.batchit(ret)
