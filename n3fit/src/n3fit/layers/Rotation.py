"""
    This module contains rotations between basis
"""

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class Rotation(MetaLayer):
    """
    Applies a transformation from the dimension-8 evolution basis
    to the dimension-14 evolution basis used by the fktables.

    The input to this layer is a `pdf_raw` variable which is expected to have
    a shape (1,  None, 8), and it is then rotated to an output (1, None, 14)
    """
    # TODO: Generate a rotation matrix in the input and just do tf.tensordot in call
    # the matrix should be: (8, 14) so that we can just do tf.tensordot(pdf, rotmat, axes=1)
    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super().__init__(**kwargs, name="evolution")

    def call(self, pdf_raw):
        # Transpose the PDF so that the flavour index is the first one
        x = op.transpose(pdf_raw)
        pdf_raw_list = [
            0 * x[0],  # photon
            x[0],  # sigma
            x[1],  # g
            x[2],  # v
            x[3],  # v3
            x[4],  # v8
            x[2],  # v15
            x[2],  # v24
            x[2],  # v35
            x[5],  # t3
            x[6],  # t8
            x[0] - 4 * x[7],  # t15 (c-)
            x[0],  # t24
            x[0],  # t35
        ]
        ret = op.concatenate(pdf_raw_list)
        # Concatenating destroys the batch index so we have to regenerate it
        return op.batchit(ret)
