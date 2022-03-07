"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.
"""

import numpy as np
from .observable import Observable
from n3fit.backends import operations as op


class DIS(Observable):
    """
        The DIS class receives a list of active flavours and a fktable
        and prepares a layer that performs the convolution of said fktable with
        the incoming pdf.

        The fktable is expected to be rank 3 (ndata, xgrid, flavours)
        while the input pdf is rank 4 where the first dimension is the batch dimension
        and the last dimension the number of replicas being fitted (1, xgrid, flavours, replicas)
    """

    def gen_mask(self, basis):
        """
            Receives a list of active flavours and generates a boolean mask tensor

            Parameters
            ----------
                basis: list(int)
                    list of active flavours
        """
        # TODO: Check the below as `basis_mask` could be unbound.
        if basis is None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def call(self, pdf):
        """
            This function perform the fktable \otimes pdf convolution.

            First pass the input PDF through a mask to remove the unactive flavors,
            then a tensor_product between the PDF and each fktable is performed
            finally the defined operation is applied to all the results

            Parameters
            ----------
                pdf:  backend tensor
                    rank 4 tensor (batch_size, xgrid, flavours, replicas)

            Returns
            -------
                result: backend tensor
                    rank 3 tensor (batchsize, replicas, ndata)
        """
        # DIS never needs splitting
        if self.splitting is not None:
            raise ValueError("DIS layer call with a dataset that needs more than one xgrid?")

        results = []
        # Separate the two possible paths this layer can take
        # TODO: Check which one is better between `op.split` and `op.gather`
        # TODO: Possibly merge these two separtations since both are now in a loop
        if self.many_masks:
            for a, (mask, fktable) in enumerate(zip(self.all_masks, self.fktables)):
                # Get the values of A & Z from the dictionary. If None, the default
                # is the free-proton ie. A=Z=1.
                a_value = 1 if self.a_list is None else self.a_list[a]
                z_value = 1 if self.z_list is None else self.z_list[a]

                # Check index of the nucleus/proton in the ordered list and select
                # only that PDF for the convolution.
                output_index = self.map_pdfs.index(a_value)
                split_size = int(pdf.shape[2] / self.nfl)
                splitted_pdf = op.split(pdf, num_or_size_splits=split_size, axis=2)
                output_pdf_range = splitted_pdf[output_index]

                # Perform standard convolution to either proton or nuclear PDF
                pdf_masked = op.boolean_mask(output_pdf_range, mask, axis=2)
                res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                results.append(res)
        else:
            for a, fktable in enumerate(self.fktables):
                # Get the values of A & Z from the dictionary. If None, the default
                # is the free-proton ie. A=Z=1.
                a_value = 1 if self.a_list is None else self.a_list[a]
                z_value = 1 if self.z_list is None else self.z_list[a]

                # Check index of the nucleus/proton in the ordered list and select
                # only that PDF for the convolution.
                output_index = self.map_pdfs.index(a_value)
                split_size = int(pdf.shape[2] / self.nfl)
                splitted_pdf = op.split(pdf, num_or_size_splits=split_size, axis=2)
                output_pdf_range = splitted_pdf[output_index]

                # Perform standard convolution to either proton or nuclear PDF
                pdf_masked = op.boolean_mask(output_pdf_range, self.all_masks[0], axis=2)
                res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                results.append(res)

        return self.operation(results)
