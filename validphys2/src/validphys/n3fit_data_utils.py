"""
n3fit_data_utils.py

Library of function that read validphys object and make them into clean numpy objects
for its usage for the fitting framework
"""
from itertools import zip_longest
import dataclasses


def _mask_fk(fktables, fk_datas, mask):
    """Mask a dataset taking into account their protections status"""
    ret = []
    for fk, fk_data in zip(fktables, fk_datas):
        if fk_data.protected:
            # TODO: this cannot be correct?
            ret.append(fk)
        else:
            ret.append(fk[mask])
    return ret


@dataclasses.dataclass
class FittableDataSet:
    """
    Python version of the libNNPDF dataset
    to be merged with the product of the new CommonData dataset
    """

    # TODO:
    # this class is basically equal to the normal DataSet
    # plus calls to generate a training_fktable, validation_fktable, etc
    # _once_ the python commondata is completed this should inherit from
    # the python-dataset adding the masking methods

    name: str
    # Now, this is a confusing one because I want to hold
    # the list of all FKTableData
    # but I want to also have methods to return directly a masked fktable
    # so I'm not sure what the names should be...
    fktables_data: list  # of validphys.coredata.FKTableData objects
    # Things that can have default values:
    operation: str = "NULL"
    frac: float = 1.0
    training_mask: bool = None

    def __post_init__(self):
        self._raw_fktables = None

    @property
    def ndata(self):
        return self.fktables_data[0].ndata

    @property
    def hadronic(self):
        return self.fktables_data[0].hadronic

    def fktables(self):
        if self._raw_fktables is None:
            self._raw_fktables = [i.get_np_fktable() for i in self.fktables_data]
        return self._raw_fktables

    def training_fktables(self):
        if self.training_mask is not None:
            return _mask_fk(self.fktables(), self.fktables_data, self.training_mask)
        return self.fktables()

    def validation_fktables(self):
        if self.training_mask is not None:
            return _mask_fk(self.fktables(), self.fktables_data, ~self.training_mask)
        return self.fktables()

    def set_mask(self, mask):
        self.training_mask = mask


def validphys_group_extractor(datasets, tr_masks):
    """
    Receives a groupping spec from validphys (most likely and experiment)
    and loops over its content extracting and parsing all information required for the fit

    Parameters
    ----------
        datasets: list of DataSetSpecs
        tr_masks: list of training mask to be applied to each dataset

    Returns
    -------
        parsed_observables: a list of (for now dictionaries) containing the information
                              required to fit the given observable
    """
    ret = []
    for dataset_spec, mask in zip_longest(datasets, tr_masks):
        # Load all fktables
        fktables = [fk.load_with_cuts(dataset_spec.cuts) for fk in dataset_spec.fkspecs]
        # And now put them in an object that contains the same information as dataset_spec save for the fact that the tables have been loaded
        # TODO: we could have a `DataSetSpec.load_for_fit(mask)` method I guess
        ret.append(
            FittableDataSet(dataset_spec.name, fktables, dataset_spec.op, dataset_spec.frac, mask)
        )
    return ret
