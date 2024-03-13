from abc import ABC, abstractmethod

import numpy as np

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

from validphys.pdfgrids import xplotting_grid


def is_unique(list_of_arrays):
    """Check whether the list of arrays more than one different arrays"""
    the_first = list_of_arrays[0]
    for i in list_of_arrays[1:]:
        if not np.array_equal(the_first, i):
            return False
    return True


def compute_pos_boundary(pdf, q0_value, xgrid, n_std, n_replicas):
    """
    Computes the Unpolarized PDF set to be used as a Boundary Condition and
    convert it into a Tensor object that can be understood by the convolution.

    Parameters
    ----------
    pdf: validphys.core.PDF
        a validphys instance of the Unpolarized PDF set
    q0_value: float
        starting scale of the theory as defined in the FK tables
    xgrid: np.ndarray
        a grid containing the x-values to be given as input to the PDF
    n_std: int
        integer representing the shift to the CV w.r.t. the standard
        deviation
    n_replicas: int
        number of replicas fitted simultaneously

    Returns
    -------
    tf.tensor:
        a tensor object that has the same shape of the output of the NN
    """
    xpdf_obj = xplotting_grid(pdf, q0_value, xgrid, basis="FK_BASIS")
    # Transpose: take the shape from (n_fl, n_x) -> (n_x, n_fl)
    xpdf_cvs = xpdf_obj.grid_values.central_value().T
    xpdf_std = xpdf_obj.grid_values.std_error().T

    # Computes the shifted Central Value as given by `n_std`
    xpdf_bound = xpdf_cvs + n_std * xpdf_std

    # Expand dimensions for multi-replicas and convert into tensor
    mult_resx = np.repeat([xpdf_bound], n_replicas, axis=0)
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
        positivity_bound,
        operation_name,
        nfl=14,
        n_replicas=1,
        **kwargs
    ):
        super(MetaLayer, self).__init__(**kwargs)

        self.dataset_name = dataset_name
        self.nfl = nfl
        self.fitbasis = fitbasis
        self.nfks = len(fktable_data)

        self.is_polarised_fktable = []  # List[bool] for polarised FK tables
        self.boundary_pdf = [None] * len(fktable_data)
        self.num_replicas = n_replicas  # TODO: Is there a reason this had to be in build?
        self.compute_observable = None  # A function (pdf, padded_fk) -> observable set in build

        all_bases = []
        xgrids = []
        fktables = []
        for idx, (fkdata, fk) in enumerate(zip(fktable_data, fktable_arr)):
            xgrids.append(fkdata.xgrid.reshape(1, -1))
            all_bases.append(fkdata.luminosity_mapping)
            fktables.append(op.numpy_to_tensor(fk))
            self.is_polarised_fktable.append(fkdata.is_polarized)

            if self.is_polarised_posdata() and fkdata.is_polarized:
                self.boundary_pdf[idx] = compute_pos_boundary(
                    pdf=positivity_bound["unpolarized_bc"],
                    q0_value=fkdata.Q0,
                    xgrid=fkdata.xgrid,
                    n_std=positivity_bound.get("n_std", 0.0),
                    n_replicas=1,  # TODO: fix this
                )
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
        for idx, (pdf, padded_fk) in enumerate(zip(pdfs, self.padded_fk_tables)):
            pdf_to_convolute = pdf if self.boundary_pdf[idx] is None else self.boundary_pdf[idx]
            observable = self.compute_observable(pdf_to_convolute, padded_fk)
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
