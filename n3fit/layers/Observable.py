import sys
from backends import MetaLayer


class Observable(MetaLayer):
    """
    This class is the parent of the DIS and DY convolutions.
    All backend-dependent code necessary for the convolutions
                                is (must be) concentrated here

    The methods gen_basis and call must be overriden by the observables
    where
        - gen_basis: it is called by the initializer and generates the mask between
                     fktables and pdfs
        - call: this is what does the actual operation
    """

    def __init__(self, output_dim, fktable, basis=None, nfl=14, **kwargs):
        self.nfl = nfl

        self.output_dim = output_dim
        self.fktable = self.np_to_tensor(fktable)
        self.xgrid_size = self.fktable.shape[-1]

        self.gen_basis(basis)

        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (self.output_dim, None)

    # Overridables
    def gen_basis(self, basis):
        print("{0} must implement the method gen_basis".format(self))
        sys.exit(-1)

    def call(self, pdf_in):
        print("{0} must implement the method gen_basis".format(self))
        sys.exit(-1)
