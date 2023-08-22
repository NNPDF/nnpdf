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
        rotation_axis: int
            rotation_axis of input to be rotated
    """

    def __init__(self, rotation_matrix, rotation_axis=3, **kwargs):
        self.rotation_matrix = op.numpy_to_tensor(rotation_matrix)
        self.rotation_axis = rotation_axis
        super().__init__(**kwargs)

    def is_identity(self):
        """Returns true if the rotation is an identity"""
        # check whether it is a mxm matrix
        if self.rotation_matrix.shape[0] == self.rotation_matrix.shape[1]:
            # check whether it is the identity
            iden = np.identity(self.rotation_matrix.shape[0])
            return np.allclose(self.rotation_matrix, iden)

    def call(self, x_raw):
        rotated = op.tensor_product(x_raw, self.rotation_matrix, [self.rotation_axis, 0])
        # this puts the rotated axis back in the original place
        return op.swapaxes(rotated, -1, self.rotation_axis)


class FlavourToEvolution(Rotation):
    """
    Rotates from the flavour basis to
    the evolution basis.
    """

    def __init__(
        self,
        flav_info,
        fitbasis,
        **kwargs,
    ):
        rotation_matrix = pdfbases.fitbasis_to_NN31IC(flav_info, fitbasis)
        super().__init__(rotation_matrix, **kwargs)


class FkRotation(Rotation):
    """
    Applies a transformation from the dimension-9 evolution basis
    to the dimension-14 evolution basis used by the fktables.

    The input to this layer is a `pdf_raw` variable which is expected to have
    a shape (1,  None, 9), and it is then rotated to an output (1, None, 14)
    """

    def __init__(self, output_dim=14, name="evolution", **kwargs):
        self.output_dim = output_dim
        rotation_matrix = self._create_rotation_matrix()
        super().__init__(rotation_matrix, name=name, **kwargs)

    def _create_rotation_matrix(self):
        """Create the rotation matrix"""
        array = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 0],  # photon
                [1, 0, 0, 0, 0, 0, 0, 0, 0],  # sigma
                [0, 1, 0, 0, 0, 0, 0, 0, 0],  # g
                [0, 0, 1, 0, 0, 0, 0, 0, 0],  # v
                [0, 0, 0, 1, 0, 0, 0, 0, 0],  # v3
                [0, 0, 0, 0, 1, 0, 0, 0, 0],  # v8
                [0, 0, 0, 0, 0, 0, 0, 0, 1],  # v15
                [0, 0, 1, 0, 0, 0, 0, 0, 0],  # v24
                [0, 0, 1, 0, 0, 0, 0, 0, 0],  # v35
                [0, 0, 0, 0, 0, 1, 0, 0, 0],  # t3
                [0, 0, 0, 0, 0, 0, 1, 0, 0],  # t8
                [1, 0, 0, 0, 0, 0, 0, -4, 0],  # t15 (c-)
                [1, 0, 0, 0, 0, 0, 0, 0, 0],  # t24
                [1, 0, 0, 0, 0, 0, 0, 0, 0],  # t35
            ]
        )
        tensor = op.numpy_to_tensor(array.T)
        return tensor


class AddPhoton(MetaLayer):
    """
    Changes the value of the photon component of the PDF to non-zero.
    The photon idx in the dimension-14 PDF basis of the FKTables is always index 0.

    In order to avoid bottlenecks, this layer can only compute the photon
    for a given fixed shape.
    In order to change the shape it is necessary to rebuild the photon.
    """

    def __init__(self, photons, **kwargs):
        self._photons_generator = photons
        self._pdf_ph = None
        super().__init__(**kwargs)

    def register_photon(self, xgrid):
        """Compute the photon array and set the layer to be rebuilt"""
        if self._photons_generator:
            self._pdf_ph = self._photons_generator(xgrid)
            self.built = False

    def call(self, pdfs, ph_replica):
        if self._pdf_ph is None:
            return pdfs
        return op.concatenate([self._pdf_ph[ph_replica], pdfs[:, :, 1:]], axis=-1)


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
