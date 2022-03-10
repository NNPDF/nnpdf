"""
n3fit_data_utils.py

This module reads validphys :py:class:`validphys.core.DataSetSpec`
and extracts the relevant information into :py:class:`validphys.n3fit_data_utils.FittableDataSet`

The ``validphys_group_extractor`` will loop over every dataset of a given group
loading their fktables (and applying any necessary cuts).
"""
from itertools import zip_longest
import dataclasses


def _mask_fk(fktables, fk_datas, mask):
    """Mask a dataset taking into account their protections status"""
    ret = []
    for fk, fk_data in zip(fktables, fk_datas):
        if fk_data.protected:
            # The protected datasets are always 1-point datasets
            # used as denominator, so if we got here this is safe
            ret.append(fk)
        else:
            ret.append(fk[mask])
    return ret


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
    # This class extends the libNNPDF dataset class
    # after the class is moved to python, this can inherit from dataset
    # and the dataset should know how to generate its "fittable" version

    name: str
    fktables_data: list  # of validphys.coredata.FKTableData objects

    # Things that can have default values:
    operation: str = "NULL"
    frac: float = 1.0
    training_mask: bool = None

    def __post_init__(self):
        self._raw_fktables = None

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
        if self._raw_fktables is None:
            self._raw_fktables = [i.get_np_fktable() for i in self.fktables_data]
        return self._raw_fktables

    def training_fktables(self):
        """Return the fktable tensors for the trainig data"""
        if self.training_mask is not None:
            return _mask_fk(self.fktables(), self.fktables_data, self.training_mask)
        return self.fktables()

    def validation_fktables(self):
        """Return the fktable tensors for the validation data"""
        if self.training_mask is not None:
            return _mask_fk(self.fktables(), self.fktables_data, ~self.training_mask)
        return self.fktables()


def validphys_group_extractor(datasets, tr_masks):
    """
    Receives a groupping spec from validphys (most likely and experiment)
    and loops over its content extracting and parsing all information required for the fit

    Parameters
    ----------
        datasets: list(:py:class:`validphys.core.DataSetSpec`)
            List of dataset spec in this group
        tr_masks: list(np.array)
            List of training mask to be set for each dataset

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
