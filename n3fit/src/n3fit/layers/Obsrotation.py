from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

class ObsRotation(MetaLayer):
    """
    Rotation is a layer used to apply a rotation transformation
    input transform matrix needs to be np array of N_out*N_in so when the
    matrix multiplication has taken place you get N_out, ... tensor out.
    If input is a true rotation then N_out=N_in
    """

    def __init__(self, transform_matrix, **kwargs):
        self.rotation = op.numpy_to_tensor(transform_matrix)
        super(MetaLayer, self).__init__(**kwargs)

    def call(self, prediction_in):
        pinT = op.transpose(prediction_in)
        return op.tensor_product(self.rotation, pinT, axes=1)
