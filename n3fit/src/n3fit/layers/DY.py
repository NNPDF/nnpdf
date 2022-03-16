import numpy as np
from .observable import Observable
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
                rank 4 tensor (batchsize, xgrid, flavours, replicas)

        Returns
        -------
            results: tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        # Hadronic observables might need splitting of the input pdf in the x dimension
        # so we have 3 different paths for this layer

        results = []
        # TODO: Check which one is better between `op.split` and `op.gather`
        if self.many_masks:
            if self.splitting:
                # TODO: Avoid repeated computations by improving the loop below
                for a, (mask, fk) in enumerate(zip(self.all_masks, self.fktables)):
                    # First we need to split the PDFs to select the corresponding f^A
                    a_value = 1 if self.a_list is None else self.a_list[a]
                    z_value = 1 if self.z_list is None else self.z_list[a]
                
                    output_index = self.map_pdfs.index(a_value)
                    split_size = int(pdf_raw.shape[2] / self.nfl)
                    splitted_pdf = op.split(pdf_raw, num_or_size_splits=split_size, axis=2)
                    nonsplitted_pdf1 = splitted_pdf[0]
                    nonsplitted_pdf2 = splitted_pdf[output_index]

                    # Then performing the standard splitting of a given f^A
                    pdf1 = op.split(nonsplitted_pdf1, self.splitting, axis=1)[a]
                    # Notice that the nuclear PDF below is only the bound-proton
                    pdf2 = op.split(nonsplitted_pdf2, self.splitting, axis=1)[a]

                    # Compute bound-neutron PDF out of the bound-proton
                    # TODO: offload some of the computations below elsewhere. The way it is done
                    # now gives myself a nightmare. At least skip them for free-proton fit.
                    if a_value != z_value:
                        neutron_pdf = op.extract_neutron_pdf(pdf2, self.neutron_mask)
                        # Compute nulcear/proton PDF out of the bound-neutron/proton PDFs
                        pdf2 = z_value * pdf2 + (a_value - z_value) * neutron_pdf
                        pdf2 /= a_value

                    pdf_x_pdf = op.pdf_masked_convolution(pdf1, pdf2, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                    results.append(res)
            else:
                for a, (mask, fk) in enumerate(zip(self.all_masks, self.fktables)):
                    # First we need to split the PDFs to select the corresponding f^A
                    a_value = 1 if self.a_list is None else self.a_list[a]
                    z_value = 1 if self.z_list is None else self.z_list[a]
                
                    output_index = self.map_pdfs.index(a_value)
                    split_size = int(pdf_raw.shape[2] / self.nfl)
                    splitted_pdf = op.split(pdf_raw, num_or_size_splits=split_size, axis=2)

                    # Then performing the standard splitting of a given f^A
                    pdf1 = splitted_pdf[0]
                    # Notice that the nuclear PDF below is only the bound-proton
                    pdf2 = splitted_pdf[output_index]

                    # Compute bound-neutron PDF out of the bound-proton
                    # TODO: offload some of the computations below elsewhere. The way it is done
                    # now gives myself a nightmare. At least skip them for free-proton fit.
                    if a_value != z_value:
                        neutron_pdf = op.extract_neutron_pdf(pdf2, self.neutron_mask)
                        # Compute nulcear/proton PDF out of the bound-neutron/proton PDFs
                        pdf2 = z_value * pdf2 + (a_value - z_value) * neutron_pdf
                        pdf2 /= a_value

                    pdf_x_pdf = op.pdf_masked_convolution(pdf1, pdf2, mask)
                    res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                    results.append(res)
        else:
            for a, fk in enumerate(self.fktables):
                # First we need to split the PDFs to select the corresponding f^A
                a_value = 1 if self.a_list is None else self.a_list[a]
                z_value = 1 if self.z_list is None else self.z_list[a]
                
                output_index = self.map_pdfs.index(a_value)
                split_size = int(pdf_raw.shape[2] / self.nfl)
                splitted_pdf = op.split(pdf_raw, num_or_size_splits=split_size, axis=2)

                # Then performing the standard splitting of a given f^A
                pdf1 = splitted_pdf[0]
                # Notice that the nuclear PDF below is only the bound-proton
                pdf2 = splitted_pdf[output_index]

                # Compute bound-neutron PDF out of the bound-proton
                # TODO: offload some of the computations below elsewhere. The way it is done
                # now gives myself a nightmare. At least skip them for free-proton fit.
                if a_value != z_value:
                    neutron_pdf = op.extract_neutron_pdf(pdf2, self.neutron_mask)
                    # Compute nulcear/proton PDF out of the bound-neutron/proton PDFs
                    pdf2 = z_value * pdf2 + (a_value - z_value) * neutron_pdf
                    pdf2 /= a_value

                pdf_x_pdf = op.pdf_masked_convolution(pdf1, pdf2, self.all_masks[0])
                res = op.tensor_product(fk, pdf_x_pdf, axes=3)
                results.append(res)

        # the masked convolution removes the batch dimension
        ret = op.transpose(self.operation(results))
        return op.batchit(ret)
