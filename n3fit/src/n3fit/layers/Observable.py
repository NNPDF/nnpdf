from n3fit.backends import MetaLayer
import tensorflow as tf
import numpy as np
from abc import abstractmethod, ABC
from n3fit.backends import operations as op

def npset(list_of_arrays):
    aa = [list_of_arrays[0]]
    for i in list_of_arrays[1:]:
        if not np.array_equal(aa[0], i):
            aa.append(i)
    return aa




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


        Parameters
        ----------
            fktabke_dicts: list(dict)
                list of fkttable_dictionaries which need to contain at least
                the 'fktable' and 'basis' entries
            operation_name: str
                string defining the name of the operation to be applied to the fktables
            nfl: int
                number of flavours in the pdf (default:14)
    """

    def __init__(self, fktable_dicts, operation_name, nfl=14, **kwargs):
        self.nfl = nfl

        basis = []
        xgrids = []
        self.fktables = []
        for fktable in fktable_dicts:
            xgrids.append(fktable['xgrid'])
            basis.append(fktable['basis'])
            self.fktables.append(op.numpy_to_tensor(fktable['fktable']))

        # Since xgrids and basis are usually 1 repeated X times let's capture
        # that behaviour
        if len(npset(xgrids)) == len(npset(basis)) == 1:
            self.basis = [self.gen_basis(basis[0])]
            self.splitting  = None
        else: # if at least one is different, they all are
            self.basis = [self.gen_basis(i) for i in basis]
            self.splitting = [i.shape[1] for i in xgrids]

        self.n_conv = len(self.fktables)
        self.operation = op.c_to_py_fun(operation_name)

        self.output_dim = self.fktables[0].shape[0]

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
    def meta_call(self, pdf_in):
        pass
