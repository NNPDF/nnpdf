"""
n3fit_data_utils.py

This module reads validphys :py:class:`validphys.core.DataSetSpec`
and extracts the relevant information into :py:class:`validphys.n3fit_data_utils.FittableDataSet`

The ``validphys_group_extractor`` will loop over every dataset of a given group
loading their fktables (and applying any necessary cuts).
"""

import dataclasses


@dataclasses.dataclass
class FittableDataSet:
    """
    Representation of the DataSet information necessary to run a fit

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

    name: str
    fktables_data: list  # of validphys.coredata.FKTableData objects

    # Things that can have default values:
    operation: str = "NULL"

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


def validphys_group_extractor(datasets):
    """
    Receives a grouping spec from validphys (most likely an experiment)
    and loops over its content extracting and parsing all information required for the fit

    Parameters
    ----------
        datasets: list(:py:class:`validphys.core.DataSetSpec`)
            List of dataset specs in this group

    Returns
    -------
        loaded_obs: list (:py:class:`validphys.n3fit_data_utils.FittableDataSet`)
    """
    loaded_obs = []
    # Use zip_longest since tr_mask can be (and it is fine) an empty list
    for dspec in datasets:
        # Load all fktables with the appropiate cuts
        fktables = [fk.load_with_cuts(dspec.cuts) for fk in dspec.fkspecs]
        # And now put them in a FittableDataSet object which
        loaded_obs.append(FittableDataSet(dspec.name, fktables, dspec.op))
    return loaded_obs
