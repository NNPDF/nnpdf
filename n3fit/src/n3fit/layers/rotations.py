"""
    This module includes rotation layers
"""
import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations as op
from validphys import pdfbases


class Rotation(MetaLayer):
    """
    Rotates the input through some user defined rotation matrix.
    Given an input matrix M_{m,n} with an input x_{m}, returns
    y_{n} = x_{m}M_{m,n}

    Parameters
    ----------
        rotation_matrix: np.array
            rotation matrix
        axes: int or list
            if given a number, contracts as many indices as given
            if given a list (of tuples) contracts indices according to op.tensor_product
    """

    def __init__(self, rotation_matrix, axes=1, **kwargs):
        self.rotation_matrix = op.numpy_to_tensor(rotation_matrix)
        self.axes = axes
        super().__init__(**kwargs)

    def is_identity(self):
        """ Returns true if the rotation is an identity """
        # check whether it is a mxm matrix
        if self.rotation_matrix.shape[0] == self.rotation_matrix.shape[1]:
            # check whether it is the identity
            iden = np.identity(self.rotation_matrix.shape[0])
            return np.allclose(self.rotation_matrix, iden)

    def call(self, x_raw):
        return op.tensor_product(x_raw, self.rotation_matrix, self.axes)


class FlavourToEvolution(Rotation):
    """
        Rotates from the flavour basis to
        the evolution basis.
    """

    def __init__(
        self, flav_info, fitbasis, **kwargs,
    ):
        rotation_matrix = pdfbases.fitbasis_to_NN31IC(flav_info, fitbasis)
        super().__init__(rotation_matrix, axes=1, **kwargs)


class FkRotation(MetaLayer):
    """
    Applies a transformation from the dimension-8 evolution basis
    to the dimension-14 evolution basis used by the fktables.

    The input to this layer is a `pdf_raw` variable which is expected to have
    a shape (1,  None, 8), and it is then rotated to an output (1, None, 14)
    """

    # TODO: Generate a rotation matrix in the input and just do tf.tensordot in call
    # the matrix should be: (8, 14) so that we can just do tf.tensordot(pdf, rotmat, axes=1)
    # i.e., create the matrix and inherit from the Rotation layer above
    def __init__(self, output_dim=14, name="evolution", **kwargs):
        self.output_dim = output_dim
        super().__init__(name=name, **kwargs)

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


class ObsRotation(MetaLayer):
    """
    Rotation is a layer used to apply a rotation transformation
    input transform matrix needs to be np array of N_out*N_in so when the
    matrix multiplication has taken place you get N_out, ... tensor out.
    If input is a true rotation then N_out=N_in
    """

    def __init__(self, transform_matrix, **kwargs):
        self.rotation = op.numpy_to_tensor(transform_matrix.T)
        super(MetaLayer, self).__init__(**kwargs)

    def call(self, prediction_in):
        return op.tensor_product(prediction_in, self.rotation, axes=1)
