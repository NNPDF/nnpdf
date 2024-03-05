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
        """
        Combine an fk table and a mask into a masked fk table to be contracted with the full PDF.

        In the case of 1 replica, this is less efficient than masking the PDF directly, so we
        leave them separate.

        Parameters
        ----------
            fk: tensor
                FK table of shape (ndata, active_flavours, x, y)
            mask: tensor
                mask of shape (flavours, flavours, active_flavours)

        Returns
        -------
            masked_fk: tensor of shape ndata, x, flavours, y, flavours) (>1 replicas case)
            (mask, fk): tuple of inputs (1 replica case)
        """
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
            self.compute_observable = compute_dy_observable_many_replica
        else:
            self.compute_observable = compute_dy_observable_one_replica


def compute_dy_observable_many_replica(pdf, masked_fk):
    """
    Contract masked fk table with two PDFs.

    Parameters
    ----------
        pdf: tensor
            pdf of shape (batch=1, replicas, xgrid, flavours)
        masked_fk: tensor
            masked fk table of shape (ndata, xgrid, flavours, xgrid, flavours)

    Returns
    -------
        tensor
            observable of shape (batch=1, replicas, ndata)
    """
    temp = op.einsum('nxfyg, bryg -> brnxf', masked_fk, pdf)
    return op.einsum('brnxf, brxf -> brn', temp, pdf)


def compute_dy_observable_one_replica(pdf, mask_and_fk):
    """
    Same operations as above but a specialized implementation that is more efficient for 1 replica,
    masking the PDF rather than the fk table.
    """
    mask, unmasked_fk = mask_and_fk
    pdf = pdf[0][0]  # yg

    mask_x_pdf = op.tensor_product(mask, pdf, axes=[(2,), (1,)])  # Ffg, yg -> Ffy
    pdf_x_pdf = op.tensor_product(mask_x_pdf, pdf, axes=[(1,), (1,)])  # Ffy, xf -> Fyx
    # nFyx, Fyx -> n
    observable = op.tensor_product(unmasked_fk, pdf_x_pdf, axes=[(1, 2, 3), (0, 1, 2)])

    return op.batchit(op.batchit(observable))  # brn
