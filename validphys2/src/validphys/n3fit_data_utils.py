"""
n3fit_data_utils.py

This module reads validphys :py:class:`validphys.core.DataSetSpec`
and extracts the relevant information into :py:class:`validphys.n3fit_data_utils.FittableDataSet`

The ``validphys_group_extractor`` will loop over every dataset of a given group
loading their fktables (and applying any necessary cuts).
"""
from itertools import zip_longest
import dataclasses
import numpy as np


@dataclasses.dataclass
class FittableDataSet:
    """
    Python version of the libNNPDF dataset
    to be merged with the product of the new CommonData dataset

    Parameters
    ----------
        name: str
            name of the dataset
        fktables_data: list(:py:class:`validphys.coredata.FKTableData`)
            list of coredata fktable objects

        operation: str
            operation to be applied to the fktables in the dataset, default "NULL"
        frac: float
            fraction of the data to enter the training set
        training_mask: bool
            training mask to apply to the fktable
    """

    # NOTE:
    # This class tries to be compatible with the libNNPDF dataset class
    # after commondata is also in python, FittableDataSet can inherit from the vp dataset
    # which knows how to generate its "fittable" version.

    name: str
    fktables_data: list  # of validphys.coredata.FKTableData objects

    # Things that can have default values:
    operation: str = "NULL"
    frac: float = 1.0
    training_mask: np.ndarray = None # boolean array

    def __post_init__(self):
        self._tr_mask = None
        self._vl_mask = None
        if self.training_mask is not None:
            data_idx = self.fktables_data[0].sigma.index.get_level_values(0).unique()
            self._tr_mask = data_idx[self.training_mask].values
            self._vl_mask = data_idx[~self.training_mask].values

    @property
    def ndata(self):
        """Number of datapoints in the dataset"""
        return self.fktables_data[0].ndata

    @property
    def hadronic(self):
        """Returns true if this is a hadronic collision dataset"""
        return self.fktables_data[0].hadronic

    def fktables(self):
        """Return the list of fktable tensors for the dataset"""
        return [fk.get_np_fktable() for fk in self.fktables_data]

    def training_fktables(self):
        """Return the fktable tensors for the trainig data"""
        if self._tr_mask is not None:
            return [fk.with_cuts(self._tr_mask).get_np_fktable() for fk in self.fktables_data]
        return self.fktables()

    def validation_fktables(self):
        """Return the fktable tensors for the validation data"""
        if self._vl_mask is not None:
            return [fk.with_cuts(self._vl_mask).get_np_fktable() for fk in self.fktables_data]
        return self.fktables()


def validphys_group_extractor(datasets, tr_masks):
    """
    Receives a grouping spec from validphys (most likely an experiment)
    and loops over its content extracting and parsing all information required for the fit

    Parameters
    ----------
        datasets: list(:py:class:`validphys.core.DataSetSpec`)
            List of dataset specs in this group
        tr_masks: list(np.array)
            List of training masks to be set for each dataset

    Returns
    -------
        loaded_obs: list (:py:class:`validphys.n3fit_data_utils.FittableDataSet`)
    """
    loaded_obs = []
    # Use zip_longest since tr_mask can be (and it is fine) an empty list
    for dspec, mask in zip_longest(datasets, tr_masks):
        # Load all fktables with the appropiate cuts
        fktables = [fk.load_with_cuts(dspec.cuts) for fk in dspec.fkspecs]
        # And now put them in a FittableDataSet object which
        loaded_obs.append(FittableDataSet(dspec.name, fktables, dspec.op, dspec.frac, mask))
    return loaded_obs
