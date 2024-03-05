"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.

    Comment on the branching based on the number of replicas:
        This is purely for performance, masking the PDF is more efficient than padding the fk table
        for one replica, and so is tensordot over einsum.

        Some timings done on snellius using tensorflow 2.15.0 and varying these 2 factors:
            | CPU\GPU | einsum | tensordot |
            | -- | -- | -- |
            | mask pdf | -  | 92 \ 65 |
            |mask fk | 330 \ 53 \  | 177 \ 53 |
            
            These timings are all for one replica.
            
            Crucially, `einsum` is a requirement of the multireplica case, while `tensordot` gives a benefit of a factor of 2x for the single replica case.
            Since this branching is required anyhow,
             by masking the PDF for 1 replica instead of padding the fktable we get an extra factor of x2
"""

import numpy as np

from n3fit.backends import operations as op

from .observable import Observable


class DIS(Observable):
    """
    The DIS class receives a list of active flavours and a fktable
    and prepares a layer that performs the convolution of said fktable with
    the incoming pdf.

    The fktable is expected to be rank 3 (ndata, xgrid, flavours)
    while the input pdf is rank 4 of shape (batch_size, replicas, xgrid, flavours)
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
                rank 1 tensor (flavours)
        """
        if basis is None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def mask_fk(self, fk, mask):
        """
        Combine an fk table and a mask into a masked fk table to be contracted with the full PDF.

        Parameters
        ----------
            fk: tensor
                FK table of shape (ndata, active_flavours, x)
            mask: tensor
                mask of shape (flavours, active_flavours)

        Returns
        -------
            masked_fk: tensor
                masked fk table of shape ndata, x, flavours)
        """
        return op.einsum('fF, nFx -> nxf', mask, fk)

    def build(self, input_shape):
        super().build(input_shape)
        if self.num_replicas > 1:
            self.compute_observable = compute_dis_observable_many_replica
        else:
            self.compute_observable = compute_dis_observable_one_replica


def compute_dis_observable_many_replica(pdf, masked_fk):
    """
    Contract masked fk table with PDF.

    Parameters
    ----------
        pdf: tensor
            pdf of shape (batch=1, replicas, xgrid, flavours)
        masked_fk: tensor
            masked fk table of shape (ndata, xgrid, flavours)

    Returns
    -------
        tensor
            observable of shape (batch=1, replicas, ndata)
    """
    return op.einsum('brxf, nxf -> brn', pdf, masked_fk)


def compute_dis_observable_one_replica(pdf, masked_fk):
    """
    Same operations as above but a specialized implementation that is more efficient for 1 replica,
    masking the PDF rather than the fk table.
    """
    # TODO: check if that is actually true
    return op.tensor_product(pdf, masked_fk, axes=[(2, 3), (1, 2)])
