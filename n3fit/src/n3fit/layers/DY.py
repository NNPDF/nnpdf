import numpy as np
from n3fit.layers.Observable import Observable
from n3fit.backends import operations as op


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
        self.basis = op.numpy_to_tensor(basis_mask, dtype=bool)

    def call(self, pdf_in):
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
        pdf_x_pdf = op.pdf_masked_convolution(pdf_in, self.basis)
        result = op.tensor_product(self.fktable, pdf_x_pdf, axes=3)
        return result
