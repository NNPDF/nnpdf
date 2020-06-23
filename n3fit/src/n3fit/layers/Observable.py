from n3fit.backends import MetaLayer
import tensorflow as tf
from abc import abstractmethod, ABC
from n3fit.backends import operations as op


class Observable(MetaLayer, ABC):
    """
        This class is the parent of the DIS and DY convolutions.
        All backend-dependent code necessary for the convolutions
                                    is (must be) concentrated here

        The methods gen_basis and call must be overriden by the observables
        where
            - gen_basis: it is called by the initializer and generates the mask between
                        fktables and pdfs
            - call: this is what does the actual operation

        # Arguments:
            - `output_dim`: output dimension of the observable (can be read from the fktable)
            - `fktable`: fktable
            - `basis`: list of active combinations of flavours
            - `nfl` : total number of flavours in the pdf (default = 14)
    """

    def __init__(self, output_dim, fktable, basis=None, nfl=14, **kwargs):
        self.nfl = nfl

        self.output_dim = output_dim
        self.fktable = op.numpy_to_tensor(fktable)
        self.xgrid_size = self.fktable.shape[-1]

        self.gen_basis(basis)

        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (self.output_dim, None)

    def digest_pdf(self, pdf):
        return tf.squeeze(pdf, axis=0)

    # Overridables
    @abstractmethod
    def gen_basis(self, basis):
        pass

    @abstractmethod
    def call(self, pdf_in):
        pass
