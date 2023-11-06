from abc import ABC, abstractmethod

import numpy as np

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


def is_unique(list_of_arrays):
    """Check whether the list of arrays more than one different arrays"""
    the_first = list_of_arrays[0]
    for i in list_of_arrays[1:]:
        if not np.array_equal(the_first, i):
            return False
    return True


def generate_neutron_mask(number_fl):
    """Generate the mask to compute the neutron-bound PDFs from the
    proton ones. Assumming `isospin asymmetry` the relation between
    the two bound PDFs is trivial. Basically, the bound-neutron PDFs
    are extracted from the proton counterpart by adding a `minus` sign
    to T3 and V3.
    
    Parameters
    ----------
    number_fl: list
        `number_fl`-PDF in evolution basis (according to FK tables
        this should always be 14)

    Returns
    -------
    generate_neutron_mask: tf.tensor
        tensor object that defines the appropriate sign to each
        PDF flavour

    """

    # TODO: Make this modular depending on the Fitting Basis
    neutron_mask = np.ones(number_fl)
    neutron_mask[4] = -1  # replace V3 sign
    neutron_mask[9] = -1  # replace T3 sign
    return op.numpy_to_tensor(neutron_mask)


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
        fktable_data: list[validphys.coredata.FKTableData]
            list of FK which define basis and xgrid for the fktables in the list
        fktable_arr: list
            list of fktables for this observable
        operation_name: str
            string defining the name of the operation to be applied to the fktables
        A_values: list
            list of A values for this observable, if the observable is a compound one
            then the this list has the same length as the FK table list
        nfl: int
            number of flavours in the pdf (default:14)
    """

    def __init__(self, fktable_data, fktable_arr, operation_name, A_values, nfl=14, **kwargs):
        super(MetaLayer, self).__init__(**kwargs)

        self.A_values = A_values
        self.nfl = nfl

        # Given that in nDIS fit the FK tables in the numnerator & Denominator
        # (for Ratios) are always the same, only one FK table is passed.
        # TODO: to re-check as this might fail in some corner cases.
        if len(A_values) != len(fktable_data) and len(fktable_data) == 1:
            fktable_arr *= len(A_values)
            fktable_data *= len(A_values)
        elif len(A_values) == len(fktable_data):
            pass
        else:
            raise ValueError("Mismatch in number of targets and FK tables.")

        basis = []
        xgrids = []
        self.fktables = []
        for fkdata, fk in zip(fktable_data, fktable_arr):
            xgrids.append(fkdata.xgrid.reshape(1, -1))
            basis.append(fkdata.luminosity_mapping)
            self.fktables.append(op.numpy_to_tensor(fk))

        # check how many xgrids this dataset needs
        if is_unique(xgrids):
            self.splitting = None
        else:
            self.splitting = [i.shape[1] for i in xgrids]

        # check how many basis this dataset needs
        if is_unique(basis) and is_unique(xgrids):
            self.all_masks = [self.gen_mask(basis[0])]
            self.many_masks = False
        else:
            self.many_masks = True
            self.all_masks = [self.gen_mask(i) for i in basis]

        self.operation = op.c_to_py_fun(operation_name)
        self.output_dim = self.fktables[0].shape[0]

        # Generate the Masks to compute neutron-bound PDFs
        self.neutron_mask = generate_neutron_mask(nfl)

    def compute_output_shape(self, input_shape):
        return (self.output_dim, None)

    # Overridables
    @abstractmethod
    def gen_mask(self, basis):
        pass

    def call(self, pdf, A_to_idx):
        A_indices = [A_to_idx[i] for i in self.A_values]
        pdfs = [pdf[:, :, i] for i in A_indices]
        return pdfs
