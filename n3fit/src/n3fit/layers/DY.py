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

    def call(self, pdf, A_to_idx):
        """
        This function perform the fktable \otimes pdf \otimes pdf convolution.

        First uses the basis of active combinations to generate a luminosity tensor
        with only some flavours active.

        The concatenate function returns a rank-3 tensor (combination_index, xgrid, xgrid)
        which can in turn be contracted with the rank-4 fktable.

        Parameters
        ----------
            pdf: list
                list of rank 4 tensors (batchsize, xgrid, flavours, replicas)
                for different A values
            A_to_idx: dict
                dictionary that maps A values to the index in the pdf

        Returns
        -------
            results: tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        # Select the correct A values from the pdf
        pdfs = super().call(pdf, A_to_idx)

        # Hadronic observables might need splitting of the input pdf in the x dimension
        # so we have 3 different paths for this layer

        # TODO: Reflect case with nuclear-bound PDFs here
        results = []
        if self.many_masks:
            # TODO: not sure what to do here
            pdf_raw = pdfs[0]
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
            for pdf, fk in zip(pdfs, self.fktables):
                pdf_x_pdf = op.pdf_masked_convolution(pdf, self.all_masks[0])
                res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                results.append(res)

        # the masked convolution removes the batch dimension
        ret = op.transpose(self.operation(results))

        return op.batchit(ret)
