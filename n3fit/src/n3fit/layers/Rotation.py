"""
    This module contains rotations between basis
"""

import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class Rotation(MetaLayer):
    """
    Applies a transformation from the dimension-8 evolution basis
    to the dimension-14 evolution basis used by the fktables.

    The input to this layer is a `pdf_raw` variable which is expected to have
    a shape (1,  None, 8), and it is then rotated to an output (1, None, 14)
    """
    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super().__init__(**kwargs, name="evolution")

    def build(self, input_shape):
        # TODO: temporary
        rotation = np.zeros((14,8))
        position_changes = [
                (1,0), (2,1), (3,2), (4,3), (5,4),
                (6,2), (7,2), (8,2), (9,5), (10,6),
                (11,0), (12,0), (13,0)
                ]
        for i,j in position_changes:
            rotation[i,j] = 1.0
        rotation[11,7] = -4.0
        self.rotation_matrix = op.numpy_to_tensor(rotation.T)
        

    def meta_call(self, pdf_raw):
        ret = op.tensor_product(pdf_raw, self.rotation_matrix, axes = 1)
        return ret

#         # Transpose the PDF so that the flavour index is the first one
#         x = op.transpose(pdf_raw)
#         pdf_raw_list = [
#             0 * x[0],  # photon
#             x[0],  # sigma 1
#             x[1],  # g 2
#             x[2],  # v 3
#             x[3],  # v3 4
#             x[4],  # v8 5
#             x[2],  # v15 6
#             x[2],  # v24 7 
#             x[2],  # v35 8
#             x[5],  # t3 9
#             x[6],  # t8 10
#             x[0] - 4 * x[7],  # t15 (c-) 11
#             x[0],  # t24 12
#             x[0],  # t35 13
#         ]
#         ret = op.concatenate(pdf_raw_list)
#         # Concatenating destroys the batch index so we have to regenerate it
#         return op.batchit(ret)
