import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class MSR_Normalization(MetaLayer):
    """
        Applies the normalisation so that the PDF output fullfills the sum rules
    """

    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super(MSR_Normalization, self).__init__(**kwargs, name="normalizer")


    def build(self, input_shape):
        # Prepare the transformation
        constants = np.ones(self.output_dim)
        # v, v8, v15, v24, v35 = 3
        vs = [3, 5, 6, 7, 8]
        for i in vs:
            constants[i] = 3.0
        transformation_add = np.zeros((self.output_dim, input_shape[-1]))
        transformation_add[2, 0] = 1.0 # we will need singlet for the MSR
        transformation_div = np.zeros((self.output_dim, input_shape[-1]))
        transformation_div[2, 1] = 1.0
        transformation_div[3, 2] = 1.0
        transformation_div[4, 3] = 1.0
        transformation_div[5, 4] = 1.0
        for i in vs[2:]:
            transformation_div[i, 2] = 1.0
        epsilon = np.ones(self.output_dim)
        epsilon[2:9] = 0.0
        # Now save them as tensorflow constants
        self.c = op.numpy_to_tensor(constants)
        self.e = op.numpy_to_tensor(epsilon)
        self.ta = op.numpy_to_tensor(transformation_add.T)
        self.td = op.numpy_to_tensor(transformation_div.T)

    def meta_call(self, xgrid):
        """
            Receives as input a tensor with the value of the MSR for each PDF
            and returns a rank-1 tensor with the normalization factor A_i of each flavour
        """
        numerator = self.c + op.tensor_product(xgrid, self.ta, axes=1)
        denominator = self.e + op.tensor_product(xgrid, self.td, axes=1)
        return numerator / denominator

# This coplicated way of doing the transformation is equivalent to 
#         x = op.flatten(xgrid)
#         pdf_sr = op.concatenate(
#             [
#                 self.one,  # photon 0 
#                 self.one,  # sigma 1
#                 (self.one - x[0]) / x[1],  # g 2
#                 self.three / x[2],  # v  3
#                 self.one / x[3],  # v3 4
#                 self.three / x[4],  # v8  5
#                 self.three / x[2],  # v15 6
#                 self.three / x[2],  # v24 7
#                 self.three / x[2],  # v35 8
#                 self.one,  # t3 9
#                 self.one,  # t8 10
#                 self.one,  # t15 (c-)  11
#                 self.one,  # t24 12
#                 self.one,  # t35 13
#             ]
#         )
#         import ipdb
#         ipdb.set_trace()
#         return pdf_sr
