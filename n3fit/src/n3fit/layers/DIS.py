"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.
"""

import numpy as np
from .observable import Observable
from n3fit.backends import operations as op


def construct_pdf(pdf, target, m_mask, n_mask):
    """Construct the PDF definition depending on the Observable.

    TODO: For all NC DIS datasets, previously, isoscalarity has
    always been enforced. That is, only the proton bound PDF `pdf`
    is returned. Check the consistency of this!!!

    """

    # Extract the values of the atomic (mass) numbers
    a_value = target["A"]
    z_value = target["Z"]

    pdf_masked = op.boolean_mask(pdf, m_mask, axis=2)

    # Compute bound-neutron PDF out of the bound-proton and construct
    # the nuclear-PDFs to be convoluted with the FK tables.
    if a_value != z_value:
        neutron_pdf = op.extract_neutron_pdf(pdf, n_mask)
        neutron_pdf_masked = op.boolean_mask(neutron_pdf, m_mask, axis=2)
        # Compute nulcear/proton PDF out of the bound-neutron/proton PDFs
        pdf_masked = z_value * pdf_masked + (a_value - z_value) * neutron_pdf_masked
        # TODO: Check the Normalization if Applicable to ALL observables
        pdf_masked /= a_value

    return pdf_masked


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

        results = []
        # Separate the two possible paths this layer can take
        if self.many_masks:
            # In the case of nuclear fit, a DIS dataset might contain different x
            if self.splitting:
                splitted_pdf = op.split(pdf, self.splitting, axis=1)
                for mask, target, pdf, fktable in zip(self.all_masks, self.target_info, splitted_pdf, self.fktables):
                    pdf_masked = construct_pdf(pdf, target, mask, self.neutron_mask)
                    res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                    results.append(res)
            else:
                for mask, fktable in zip(self.all_masks, self.fktables):
                    pdf_masked = op.boolean_mask(pdf, mask, axis=2)
                    res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                    results.append(res)
        else:
            pdf_masked = op.boolean_mask(pdf, self.all_masks[0], axis=2)
            for fktable in self.fktables:
                res = op.tensor_product(pdf_masked, fktable, axes=[(1, 2), (2, 1)])
                results.append(res)

        return self.operation(results)
