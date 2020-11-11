from n3fit.backends import MetaLayer
import numpy as np
from abc import abstractmethod, ABC
from n3fit.backends import operations as op


def _is_unique(list_of_arrays):
    """ Check whether the list of arrays more than one different arrays """
    the_first = list_of_arrays[0]
    for i in list_of_arrays[1:]:
        if not np.array_equal(the_first, i):
            return False
    return True


class Observable(MetaLayer, ABC):
    """
        This class is the parent of the DIS and DY convolutions.
        All backend-dependent code necessary for the convolutions
                                    is (must be) concentrated here

        The methods gen_mask and call must be overriden by the observables
        where
            - gen_mask: it is called by the initializer and generates the mask between
                        fktables and pdfs
            - call: this is what does the actual operation


        Parameters
        ----------
            fktable_dicts: list
                list of fktable_dicts which define basis and xgrid for the fktables in the list
            fktable_arr: list
                list of fktables for this observable
            operation_name: str
                string defining the name of the operation to be applied to the fktables
            nfl: int
                number of flavours in the pdf (default:14)
    """

    def __init__(self, fktable_dicts, fktable_arr, operation_name, nfl=14, **kwargs):
        super(MetaLayer, self).__init__(**kwargs)

        self.nfl = nfl

        basis = []
        xgrids = []
        self.fktables = []
        for fktable, fk in zip(fktable_dicts, fktable_arr):
            xgrids.append(fktable["xgrid"])
            basis.append(fktable["basis"])
            self.fktables.append(op.numpy_to_tensor(fk))

        # check how many xgrids this dataset needs
        if _is_unique(xgrids):
            self.splitting = None
        else:
            self.splitting = [i.shape[1] for i in xgrids]

        # check how many basis this dataset needs
        if _is_unique(basis) and _is_unique(xgrids):
            self.all_masks = [self.gen_mask(basis[0])]
            self.many_masks = False
        else:
            self.many_masks = True
            self.all_masks = [self.gen_mask(i) for i in basis]

        self.operation = op.c_to_py_fun(operation_name)
        self.output_dim = self.fktables[0].shape[0]

    def compute_output_shape(self, input_shape):
        return (self.output_dim, None)

    # Overridables
    @abstractmethod
    def gen_mask(self, basis):
        pass
