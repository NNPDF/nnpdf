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


def compute_posbc(extern_pdf, n_replicas, idx):
    """
    Extract the relevant PDF to be used a Boundary Condition and convert
    it into a Tensor that can be understood by the convolution.

    extern_lhapdf: list[np.ndarray]
        list of pre-computed PDF for a fixed Q2 with shape (n_x, n_fl),
        the length of the list correspond to the number of FK tables
        required for the Positivity dataset
    n_replicas: int
        number of replicas
    idx: int
        index specifying which element of `extern_lhapdf` should be used
    """
    mult_resx = np.repeat([extern_pdf[idx]], n_replicas, axis=0)
    add_batch_resx = np.expand_dims(mult_resx, axis=0)
    return op.numpy_to_tensor(add_batch_resx)


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
        nfl=14,
        **kwargs
    ):
        super(MetaLayer, self).__init__(**kwargs)

        self.dataset_name = dataset_name
        self.nfl = nfl
        self.fitbasis = fitbasis
        self.nfks = len(fktable_data)

        self.is_polarised_fktable = [] # List[bool] for polarised FK tables
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
            self.is_polarised_fktable.append(fkdata.is_polarized)

            if self.is_polarised_posdata():
                self.pdfbc.append(extern_lhapdf(fkdata.xgrid.tolist()))
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
        # NOTE: Here, effectively, we are computing: (q + qbar) - |Δq| + |Δqbar|
        # which should be fine because: |Δq+Δqbar| <= |Δq| + |Δqbar|
        for idx, (pdf, padded_fk) in enumerate(zip(pdfs, self.padded_fk_tables)):
            # Check if Unpolarized POS FK and convolute with the pre-computed PDF
            if self.is_polarised_posdata() and not self.is_polarised_fktable[idx]:
                pdf_to_convolute = compute_posbc(self.pdfbc, self.num_replicas, idx)
            else: # Otherwsie, convolute with the usual NN PDFs
                pdf_to_convolute = pdf

            # Compute the usual convolution between the (pre-computed) PDF
            observable = self.compute_observable(pdf_to_convolute, padded_fk)

            # If it is instead a Polarized POS FK, we need to take its Absolute
            if self.is_polarised_posdata() and self.is_polarised_fktable[idx]:
                observable = op.multiply_minusone(op.absolute(observable))

            observables.append(observable)

        observables = self.operation(observables)
        return observables

    def is_polarised_posdata(self):
        """Check if the dataset is a Polarised Positivity dataset."""
        # TODO: Better move the following to the new commondata input
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
