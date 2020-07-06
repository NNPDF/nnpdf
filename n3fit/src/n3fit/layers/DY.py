import numpy as np
from n3fit.layers.Observable import Observable
from n3fit.backends import operations as op


class DY(Observable):
    """
    Computes the convolution of two PDFs (the same one twice) and one fktable
    """

    def gen_mask(self, basis):
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
                rank 2 tensor (batchsize, ndata)
        """
<<<<<<< HEAD
        pdf_in = self.digest_pdf(pdf_in)
        # This is a convoluted way of applying a mask, but it is faster
        # mask-version below
        lumi_fun = []
        pdfT = op.transpose(pdf_in)

        for i, j in self.basis:
            lumi_fun.append(op.tensor_product(pdfT[i], pdfT[j], axes=0))

        pdf_X_pdf = op.concatenate(lumi_fun, axis=0, target_shape=(self.basis_size, self.xgrid_size, self.xgrid_size))

        result = op.tensor_product(self.fktable, pdf_X_pdf, axes=3)
        return result


# Another example on how to performt the DY convolution
# this code is equivalent to the previos one, with a slightly greater cost
# class DY_mask(Observable):
#     def gen_basis(self, basis):
#         if basis is not None:
#             self.basis = np.zeros((self.nfl, self.nfl), dtype=bool)
#             for i, j in basis.reshape(-1, 2):
#                 self.basis[i, j] = True
#         else:
#             self.basis = np.ones((self.nfl, self.nfl), dtype=bool)
# 
#     def call(self, pdf_in):
#         lfun = op.tensor_product(pdf_in, pdf_in, axes=0)
#         lfunT = tensorflow.keras.backend.permute_dimensions(lfun, (3, 1, 2, 0))
#         x = op.boolean_mask(lfunT, self.basis, axis=0)
#         result = op.tensor_product(self.fktable, x, axes=3)
#         return result
=======
        # Hadronic observables might need splitting of the input pdf in the x dimension
        # so we have 3 different paths for this layer

        results = []
        if self.many_masks:
            if self.splitting:
                splitted_pdf = op.split(pdf_raw, self.splitting, axis=1)
                for mask, pdf, fk in zip(self.all_masks, splitted_pdf, self.fktables):
                    pdf_x_pdf = op.pdf_masked_convolution(pdf, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                    results.append(res)
            else:
                for mask, fk in zip(self.all_masks, self.fktables):
                    pdf_x_pdf = op.pdf_masked_convolution(pdf_raw, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                    results.append(res)
        else:
            pdf_x_pdf = op.pdf_masked_convolution(pdf_raw, self.all_masks[0])
            for fk in self.fktables:
                res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                results.append(res)

        # the masked convolution removes the batch dimension
        ret = self.operation(results)
        return op.batchit(ret)
>>>>>>> master
