import numpy as np
from n3fit.layers.Observable import Observable


class DY(Observable):
    """
    Computes the convolution of two PDFs (the same one twice) and one fktable
    """
    def gen_basis(self, basis):
        """
        Receives an array of combinations and make it into an array of 2-tuples

        # Arguments:
                - `basis`: array of combinations in the form
                           [i1,j1,i2,j2,i3,j3...]
        """
        basis_to_pairs = basis.reshape(-1, 2)
        self.basis = basis_to_pairs
        self.basis_size = len(self.basis)

    def call(self, pdf_in):
        """
            Thiss function perform the fktable \otimes pdf \otimes pdf convolution.

            First uses the basis of active combinations to generate a luminosity tensor
            with only some flavours active.

            The concatenate function returns a rank-3 tensor (combination_index, xgrid, xgrid)
            which can in turn be contracted with the rank-4 fktable.

            # Arguments:
                - `pdf_in`: rank 2 tensor (xgrid, flavours)

            # Returns:
                - `result`: rank 1 tensor (ndata)
        """
        # This is a convoluted way of applying a mask, but it is faster
        # mask-version below
        lumi_fun = []
        pdfT = self.transpose(pdf_in)

        for i, j in self.basis:
            lumi_fun.append(self.tensor_product(pdfT[i], pdfT[j], axes=0))

        pdf_X_pdf = self.concatenate(lumi_fun, axis=0, target_shape=(self.basis_size, self.xgrid_size, self.xgrid_size))

        result = self.tensor_product(self.fktable, pdf_X_pdf, axes=3)
        return result


# Another example on how to performt the DY convolution
# this code is equivalent to the previos one, with a slightly greater cost
class DY_mask(Observable):
    def gen_basis(self, basis):
        if basis is not None:
            self.basis = np.zeros((self.nfl, self.nfl), dtype=bool)
            for i, j in basis.reshape(-1, 2):
                self.basis[i, j] = True
        else:
            self.basis = np.ones((self.nfl, self.nfl), dtype=bool)

    def call(self, pdf_in):
        lfun = self.tensor_product(pdf_in, pdf_in, axes=0)
        lfunT = self.permute_dimensions(lfun, (3, 1, 2, 0))
        x = self.boolean_mask(lfunT, self.basis, axis=0)
        result = self.tensor_product(self.fktable, x, axes=3)
        return result
