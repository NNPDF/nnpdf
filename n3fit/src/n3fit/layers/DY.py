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

    def mask_fk(self, fk, mask):
        if self.num_replicas > 1:
            return op.einsum('fgF, nFxy -> nxfyg', mask, fk)
        else:
            mask = op.einsum('fgF -> Ffg', mask)
            fk = op.einsum('nFxy -> nFyx', fk)
            mask_and_fk = (mask, fk)
            return mask_and_fk

    def build(self, input_shape):
        super().build(input_shape)
        if self.num_replicas > 1:

            def compute_observable(pdf, masked_fk):
                temp = op.einsum('nxfyg, bryg -> brnxf', masked_fk, pdf)
                return op.einsum('brnxf, brxf -> brn', temp, pdf)

        else:

            def compute_observable(pdf, mask_and_fk):
                # with 1 replica, it's more efficient to mask the PDF rather than the fk table
                mask, unmasked_fk = mask_and_fk
                pdf = pdf[0][0]  # yg

                mask_x_pdf = op.tensor_product(mask, pdf, axes=[(2,), (1,)])  # Ffg, yg -> Ffy
                pdf_x_pdf = op.tensor_product(mask_x_pdf, pdf, axes=[(1,), (1,)])  # Ffy, xf -> Fyx
                # nFyx, Fyx -> n
                observable = op.tensor_product(unmasked_fk, pdf_x_pdf, axes=[(1, 2, 3), (0, 1, 2)])

                return op.batchit(op.batchit(observable))  # brn

        self.compute_observable = compute_observable
