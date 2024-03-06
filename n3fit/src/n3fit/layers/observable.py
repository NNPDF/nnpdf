from abc import ABC, abstractmethod

import numpy as np

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


# The following maps the indices to either the Unpolarized or
# Polarized Positivity FK tables.
POS_POLSD_INDEX = [0, 2]  # Polarised POS FK tables
POS_UNPOL_INDEX = [1, 3]  # Unpolarized POS FK tables


def is_unique(list_of_arrays):
    """Check whether the list of arrays more than one different arrays"""
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
        fktable_data: list[validphys.coredata.FKTableData]
            list of FK which define basis and xgrid for the fktables in the list
        fktable_arr: list
            list of fktables for this observable
        operation_name: str
            string defining the name of the operation to be applied to the fktables
        nfl: int
            number of flavours in the pdf (default:14)
    """

    def __init__(
        self,
        fktable_data,
        fktable_arr,
        dataset_name,
        fitbasis,
        extern_lhapdf,
        operation_name,
        n_replicas,
        nfl=14,
        **kwargs
    ):
        super(MetaLayer, self).__init__(**kwargs)

        self.dataset_name = dataset_name
        self.nfl = nfl
        self.fitbasis = fitbasis
        self.nfks = len(fktable_data)

        self.pdfbc = [] # Pre-computed Unpolarized PDF Boundary Condition
        self.num_replicas = None  # set in build
        self.compute_observable = None  # A function (pdf, padded_fk) -> observable set in build

        all_bases = []
        xgrids = []
        fktables = []
        for fkdata, fk in zip(fktable_data, fktable_arr):
            xgrids.append(fkdata.xgrid.reshape(1, -1))
            all_bases.append(fkdata.luminosity_mapping)
            fktables.append(op.numpy_to_tensor(fk))

            if self.is_polarised_pos():
                resx = extern_lhapdf(fkdata.xgrid.tolist())
                mult_resx = np.repeat([resx], n_replicas, axis=0)
                resx = np.expand_dims(mult_resx, axis=0)
                self.pdfbc.append(op.numpy_to_tensor(resx))
        self.fktables = fktables

        # check how many xgrids this dataset needs
        if is_unique(xgrids):
            self.splitting = None
        else:
            self.splitting = [i.shape[1] for i in xgrids]

        self.operation = op.c_to_py_fun(operation_name)
        self.output_dim = self.fktables[0].shape[0]

        if is_unique(all_bases) and is_unique(xgrids):
            self.all_masks = [self.gen_mask(all_bases[0])]
        else:
            self.all_masks = [self.gen_mask(basis) for basis in all_bases]

        self.masks = [compute_float_mask(bool_mask) for bool_mask in self.all_masks]

    def build(self, input_shape):
        self.num_replicas = input_shape[1]

        # repeat the masks if necessary for fktables (if not, the extra copies
        # will get lost in the zip)
        masks = self.masks * len(self.fktables)
        self.padded_fk_tables = [self.pad_fk(fk, mask) for fk, mask in zip(self.fktables, masks)]

        super().build(input_shape)

    def call(self, pdf):
        """
        This function perform the convolution with the fktable and one (DIS) or two (DY-like) pdfs.

        Parameters
        ----------
            pdf:  backend tensor
                rank 4 tensor (batch_size, replicas, xgrid, flavours)

        Returns
        -------
            observables: backend tensor
                rank 3 tensor (batchsize, replicas, ndata)
        """
        if self.splitting:
            pdfs = op.split(pdf, self.splitting, axis=2)
        else:
            pdfs = [pdf] * len(self.padded_fk_tables)

        observables = []
        # In the case of Polarized fits, the computation of Positivity observables
        # are much more involved. In the case of (anti-)quarks Positivity, we need
        # to compute as final observable:
        #      O = (q + qbar) - |Δq+Δqbar|   where Δq=Polarized quark PDF
        # And in the case of Gluon Positivity, we need to compute:
        #      O = g - |Δg|
        # This is why the POS FK tables are ordered in an alternating way between
        # Polarized and Unpolarized in the Metdata and hence the ordering of both
        # `POS_POLSD_INDEX` and `POS_UNPOL_INDEX`.
        for idx, (pdf, padded_fk) in enumerate(zip(pdfs, self.padded_fk_tables)):
            # Check if Unpolarized POS FK and convolute with the pre-computed PDF
            if self.is_polarised_pos() and idx in POS_UNPOL_INDEX:
                pdf_to_convolute = self.pdfbc[idx]
            else: # Convolute with the NN PDFs
                pdf_to_convolute = pdf

            # Compute the usual convolution between the (pre-computed) PDF
            observable = self.compute_observable(pdf_to_convolute, padded_fk)

            # If it is instead an Unpolarized POS FK, we need to take its Absolute
            # NOTE: Here, effectively, we are computing: (q + qbar) - |Δq| + |Δqbar|
            # which should be fine because: |Δq+Δqbar| <= |Δq| + |Δqbar|
            if self.is_polarised_pos() and idx in POS_POLSD_INDEX:
                observable = op.multiply_minusone(op.absolute(observable))

            observables.append(observable)

        observables = self.operation(observables)
        return observables

    def is_polarised_pos(self):
        if "POL" in self.fitbasis and "_POS_" in self.dataset_name:
            # Polarised POS contains at least 2 FK tables
            return self.nfks >= 2
        return False

    # Overridables
    @abstractmethod
    def gen_mask(self, basis):
        pass

    @abstractmethod
    def pad_fk(self, fk, mask):
        pass


def compute_float_mask(bool_mask):
    """
    Compute a float form of the given boolean mask, that can be contracted over the full flavor
    axes to obtain a PDF of only the active flavors.

    Parameters
    ----------
        bool_mask: boolean tensor
            mask of the active flavours

    Returns
    -------
        masked_to_full: float tensor
            float form of mask
    """
    # Create a tensor with the shape (**bool_mask.shape, num_active_flavours)
    masked_to_full = []
    for idx in np.argwhere(bool_mask):
        temp_matrix = np.zeros(bool_mask.shape)
        temp_matrix[tuple(idx)] = 1
        masked_to_full.append(temp_matrix)
    masked_to_full = np.stack(masked_to_full, axis=-1)
    masked_to_full = op.numpy_to_tensor(masked_to_full)

    return masked_to_full
