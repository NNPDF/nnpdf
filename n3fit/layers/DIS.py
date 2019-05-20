import numpy as np
from layers.Observable import Observable

class DIS(Observable):
    def gen_basis(self, basis):
        if basis is not None:
            self.basis = np.zeros(self.nfl, dtype = bool)
            for i in basis:
                self.basis[i] = True
        else:
            self.basis = np.ones(self.nfl, dtype = bool)

    def call(self, pdf_in):

        pdf_masked = self.boolean_mask(pdf_in, self.basis, axis = 1)

        pdfT = self.transpose(pdf_masked)

        result = self.tensor_product(self.fktable, pdfT, axes = 2)
        return result
