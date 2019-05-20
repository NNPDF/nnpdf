import numpy as np
from layers.Observable import Observable


class DY(Observable):
    def gen_basis(self, basis):
        basis_to_pairs = basis.reshape(-1,2)
        self.basis = basis_to_pairs
        self.basis_size = len(self.basis)

    def call(self, pdf_in):
        # This is a convoluted way of applying a mask, but it is faster
        # mask-version below
        lumi_fun = []
        pdfT = self.transpose(pdf_in)

        for i,j in self.basis:
            lumi_fun.append( self.tensor_product(pdfT[i], pdfT[j], axes = 0) )

        pdfXpdf = self.concatenate( lumi_fun, axis = 0,
                target_shape = (self.basis_size, self.xgrid_size, self.xgrid_size) )

        result = self.tensor_product(self.fktable, pdfXpdf, axes = 3)
        return result

class DY_mask(Observable):
    def gen_basis(self, basis):
        if basis is not None:
            self.basis = np.zeros((self.nfl, self.nfl), dtype = bool)
            for i,j in basis.reshape(-1,2):
                self.basis[i,j] = True
        else:
            self.basis = np.ones((self.nfl, self.nfl), dtype = bool)

    def call(self, pdf_in):
        lfun = self.tensor_product(pdf_in, pdf_in, axes = 0)
        lfunT = self.permute_dimensions(lfun, (3,1,2,0))
        x = self.boolean_mask(lfunT, self.basis, axis = 0)
        result = self.tensor_product(self.fktable, x, axes = 3)
        return result
