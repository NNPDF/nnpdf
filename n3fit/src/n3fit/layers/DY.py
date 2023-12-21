import numpy as np

from n3fit.backends import operations as op

from .observable import Observable


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

    def call(self, pdf_raw):
        """
        This function perform the fktable \otimes pdf \otimes pdf convolution.

        First uses the basis of active combinations to generate a luminosity tensor
        with only some flavours active.

        The concatenate function returns a rank-3 tensor (combination_index, xgrid, xgrid)
        which can in turn be contracted with the rank-4 fktable.

        Parameters
        ----------
            pdf_in: tensor
                rank 4 tensor (batchsize, replicas, xgrid, flavours)

        Returns
        -------
            results: tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        # Hadronic observables might need splitting of the input pdf in the x dimension
        # so we have 3 different paths for this layer

        results = []
        if self.many_masks:
            if self.splitting:
                splitted_pdf = op.split(pdf_raw, self.splitting, axis=2)
                for mask, pdf, fk in zip(self.all_masks, splitted_pdf, self.fktables):
                    pdf_x_pdf = op.pdf_masked_convolution(pdf, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=[(1, 2, 3), (1, 2, 3)])
                    results.append(res)
            else:
                for mask, fk in zip(self.all_masks, self.fktables):
                    pdf_x_pdf = op.pdf_masked_convolution(pdf_raw, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=[(1, 2, 3), (1, 2, 3)])
                    results.append(res)
        else:
            pdf_x_pdf = op.pdf_masked_convolution(pdf_raw, self.all_masks[0])
            for fk in self.fktables:
                res = op.tensor_product(fk, pdf_x_pdf, axes=[(1, 2, 3), (1, 2, 3)])
                results.append(res)

        # the masked convolution removes the batch dimension
        ret = op.transpose(self.operation(results))
        return op.batchit(ret)
