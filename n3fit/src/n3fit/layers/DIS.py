import numpy as np
from n3fit.layers.Observable import Observable


class DIS(Observable):
    """
        The DIS class receives a list of active flavours and a fktable
        and prepares a layer that performs the convolution of said fktable with
        the incoming pdf.
    """
    def gen_basis(self, basis):
        """
            Receives a list of active flavours and generates a boolean mask

            # Arguments:
                - `basis`: list of active flavours
        """
        if basis is not None:
            self.basis = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                self.basis[i] = True
        else:
            self.basis = np.ones(self.nfl, dtype=bool)

    def call(self, pdf_in):
        """
            Thiss function perform the fktable \otimes pdf convolution.

            Firs pass the input PDF through a mask to remove the unactive flavour, then transpose the PDF
            to have everything in the correct order and finally perform a tensorproduct contracting both pdf indices.

            # Arguments:
                - `pdf_in`: rank 2 tensor (xgrid, flavours)
            # Returns:
                - `result`: rank 1 tensor (ndata)
        """
        pdf_masked = self.boolean_mask(pdf_in, self.basis, axis=1)

        pdfT = self.transpose(pdf_masked)
        result = self.tensor_product(self.fktable, pdfT, axes=2)
        return result
