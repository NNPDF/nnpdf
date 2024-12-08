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

    def pad_fk(self, fk, mask):
        """
        Combine an fk table and a mask into an fk table padded with zeroes for the inactive
        flavours, to be contracted with the full PDF.

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
            padded_fk: tensor of shape ndata, x, flavours, y, flavours) (>1 replicas case)
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


def compute_dy_observable_many_replica(pdf, padded_fk):
    """
    Contract masked fk table with two PDFs.

    Parameters
    ----------
        pdf: list[tensor]
            list of pdf of shape (batch=1, replicas, xgrid, flavours)
        padded_fk: tensor
            masked fk table of shape (ndata, xgrid, flavours, xgrid, flavours)

    Returns
    -------
        tensor
            observable of shape (batch=1, replicas, ndata)
    """
    pdfa = pdf[1]
    pdfb = pdf[0]

    temp = op.einsum('nxfyg, bryg -> brnxf', padded_fk, pdfa)
    return op.einsum('brnxf, brxf -> brn', temp, pdfb)


def compute_dy_observable_one_replica(pdf, mask_and_fk):
    """
    Same operations as above but a specialized implementation that is more efficient for 1 replica,
    masking the PDF rather than the fk table.
    """
    # mask: (channels, flavs_b, flavs_a) Ffg
    # fk: (npoints, channels, x_a, x_b) nFyx
    mask, fk = mask_and_fk
    # Retrieve the two PDFs (which may potentially be coming from different initial states)
    # Since this is the one-replica function, remove the batch and replica dimension
    pdfb = pdf[0][0][0]  # (x_b, flavs_b) xf
    pdfa = pdf[1][0][0]  # (x_a, flavs_a) yg

    # TODO: check which PDF must go first in case of different initial states!!!
    mask_x_pdf = op.tensor_product(mask, pdfa, axes=[(2,), (1,)])  # Ffg, yg -> Ffy
    pdf_x_pdf = op.tensor_product(mask_x_pdf, pdfb, axes=[(1,), (1,)])  # Ffy, xf -> Fyx
    observable = op.tensor_product(fk, pdf_x_pdf, axes=[(1, 2, 3), (0, 1, 2)])  # nFyx, Fyx -> n

    return op.batchit(op.batchit(observable))  # brn
