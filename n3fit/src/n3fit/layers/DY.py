import numpy as np

from n3fit.backends import operations as op

from .observable import Observable


class DY(Observable):
    """
    Computes the convolution of two PDFs (the same one twice) and one fktable
    """

    def gen_mask(self, basis):
        """
        Receives a list of active flavours and generates a boolean mask tensor

        Parameters
        ----------
            basis: list(int)
                list of active flavours

        Returns
        -------
            mask: tensor
                rank 2 tensor (flavours, flavours)
        """
        if basis is None:
            basis_mask = np.ones((self.nfl, self.nfl), dtype=bool)
        else:
            basis_mask = np.zeros((self.nfl, self.nfl), dtype=bool)
            for i, j in basis.reshape(-1, 2):
                basis_mask[i, j] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def build(self, input_shape):
        super().build(input_shape)
        if self.num_replicas > 1:

            def compute_observable(pdf, mask, fk):
                return op.einsum('fgF, nFxy, brxf, bryg -> brn', mask, fk, pdf, pdf)

        else:

            def compute_observable(pdf, mask, fk):
                pdf = pdf[0][0]  # yg
                pdf_x_mask = op.tensor_product(pdf, mask, axes=[[1], [1]])  # yg, fgF -> yfF
                pdf_x_pdf = op.tensor_product(pdf, pdf_x_mask, axes=[[1], [1]])  # xf, yfF -> xyF
                observable = op.tensor_product(
                    fk, pdf_x_pdf, axes=[(1, 2, 3), (2, 0, 1)]
                )  # nFxy, xyF
                return op.batchit(op.batchit(observable))  # brn

        self.compute_observable = compute_observable
