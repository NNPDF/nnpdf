"""
    The all-important CommonData object
"""

import dataclasses
from typing import Optional

import numpy as np
import pandas as pd

from .utils import yaml

KIN_NAMES = ["kin1", "kin2", "kin3"]


def generate_path_filtered_data(fit_path, setname):
    """Utility to ensure that both the loader and tools like setupfit utilize the same convention
    to generate the names of generated pseudodata"""
    data_path = fit_path / "filter" / setname / f"filtered_data_{setname}.yaml"
    unc_path = data_path.with_name(f"filtered_uncertainties_{setname}.yaml")
    return data_path, unc_path


@dataclasses.dataclass(eq=False)
class CommonData:
    """
    Data contained in Commondata files, relevant cuts applied.

    Parameters
    ----------

    setname : str
        Name of the dataset

    ndata : int
        Number of data points

    commondataproc : str
        Process type, one of 21 options

    nkin : int
        Number of kinematics specified

    nsys : int
        Number of systematics

    commondata_table : pd.DataFrame
        Pandas dataframe containing the commondata

    systype_table : pd.DataFrame
        Pandas dataframe containing the systype index
        for each systematic alongside the uncertainty
        type (ADD/MULT/RAND) and name
        (CORR/UNCORR/THEORYCORR/SKIP)

    systematics_table: pd.DataFrame
        Panda dataframe containing the table of systematics
    """

    setname: str
    ndata: int
    commondataproc: str
    nkin: int
    nsys: int
    commondata_table: pd.DataFrame = dataclasses.field(repr=False)
    systype_table: pd.DataFrame = dataclasses.field(repr=False)
    legacy: bool = False
    systematics_table: Optional[pd.DataFrame] = dataclasses.field(init=None, repr=False)
    legacy_names: Optional[list] = None
    kin_variables: Optional[list] = None

    def __post_init__(self):
        self.systematics_table = self.commondata_table.drop(
            columns=["process", "data", "stat"] + KIN_NAMES
        )
        # TODO: set for now commondataproc as a string
        self.commondataproc = str(self.commondataproc)

    def with_cuts(self, cuts):
        """A method to return a CommonData object where
        an integer mask has been applied, keeping only data
        points which pass cuts.

        Note if the first data point passes cuts, the first entry
        of ``cuts`` should be ``0``.

        Paramters
        ---------
        cuts: list or validphys.core.Cuts or None
        """
        # Ensure that the cuts we're applying applies to this dataset
        # only check, however, if the cuts is of type :py:class:`validphys.core.Cuts`
        if hasattr(cuts, "name") and self.setname != cuts.name:
            raise ValueError(
                f"The cuts provided are for {cuts.name} which does not apply "
                f"to this commondata file: {self.setname}"
            )

        if hasattr(cuts, "load"):
            cuts = cuts.load()
        if cuts is None:
            return self

        # We must shift the cuts up by 1 since a cut of 0 implies the first data point
        # while commondata indexing starts at 1.
        cuts = list(map(lambda x: x + 1, cuts))

        newndata = len(cuts)
        new_commondata_table = self.commondata_table.loc[cuts]
        return dataclasses.replace(self, ndata=newndata, commondata_table=new_commondata_table)

    @property
    def kinematics(self):
        return self.commondata_table[KIN_NAMES]

    def get_kintable(self):
        return self.kinematics.values

    @property
    def central_values(self):
        return self.commondata_table["data"]

    def with_central_value(self, cv):
        tb = self.commondata_table.copy()
        tb["data"] = cv
        return dataclasses.replace(self, commondata_table=tb)

    def get_cv(self):
        return self.central_values.values

    @property
    def stat_errors(self):
        return self.commondata_table["stat"]

    @property
    def multiplicative_errors(self):
        """Returns the systematics which are multiplicative (systype is MULT)
        in a percentage format, with SKIP uncertainties removed.

        """
        mult_systype = self.systype_table[self.systype_table["treatment"] == "MULT"]
        mult_table = self.systematics_table.filter(like="MULT")

        if self.legacy:
            # Needed in legacy because every uncertainty appears as both mult and add
            # so it is necessary to select the uncertainties that are to be consireded as MULT/ADD
            # Minus 1 because iloc starts from 0, while the systype counting starts from 1
            mult_table = mult_table.iloc[:, mult_systype.index - 1]

        mult_table.columns = mult_systype["name"].to_numpy()
        return mult_table.loc[:, mult_table.columns != "SKIP"]

    @property
    def additive_errors(self):
        """Returns the systematics which are additive (systype is ADD) as
        absolute uncertainties (same units as data), with SKIP uncertainties
        removed.

        """
        add_systype = self.systype_table[self.systype_table["treatment"] == "ADD"]
        add_table = self.systematics_table.filter(like="ADD")

        if self.legacy:
            # Minus 1 because iloc starts from 0, while the systype counting starts from 1
            add_table = add_table.iloc[:, add_systype.index - 1]

        add_table.columns = add_systype["name"].to_numpy()
        return add_table.loc[:, add_table.columns != "SKIP"]

    def systematic_errors(self, central_values=None):
        """Returns all systematic errors as absolute uncertainties, with a
        single column for each uncertainty. Converts
        :py:attr:`multiplicative_errors` to units of data and then appends
        onto :py:attr:`additive_errors`. By default uses the experimental
        central values to perform conversion, but the user can supply a
        1-D array of central values, with length :py:attr:`self.ndata`, to use
        instead of the experimental central values to calculate the absolute
        contribution of the multiplicative systematics.

        Parameters
        ----------
        central_values: None, np.array
            1-D array containing alternative central values to combine with
            multiplicative uncertainties. This array must have length equal
            to :py:attr:`self.ndata`. By default ``central_values`` is None, and
            the central values of the commondata are used.

        Returns
        -------
        systematic_errors: pd.DataFrame
            Dataframe containing systematic errors.

        """
        if central_values is None:
            central_values = self.central_values.to_numpy()
        converted_mult_errors = self.multiplicative_errors * central_values[:, np.newaxis] / 100
        return pd.concat((self.additive_errors, converted_mult_errors), axis=1)

    def export_data(self, buffer):
        """Exports the central data defined by this commondata instance to the given buffer"""
        ret = {"data_central": self.central_values.tolist()}
        yaml.safe_dump(ret, buffer)

    def export_uncertainties(self, buffer):
        """Exports the uncertainties defined by this commondata instance to the given buffer"""
        definitions = {}
        for idx, row in self.systype_table.iterrows():
            if row["name"] != "SKIP":
                definitions[f"sys_{idx}"] = {"treatment": row["treatment"], "type": row["name"]}

        # Order the definitions by treatment as ADD, MULT
        # TODO: make it so that it corresponds to the original order exactly
        sorted_definitions = {
            k: v for k, v in sorted(definitions.items(), key=lambda item: item[1]["treatment"])
        }
        bins = []
        for idx, row in self.systematic_errors().iterrows():
            tmp = {"stat": float(self.stat_errors[idx])}
            # Hope things come in the right order...
            for key_name, val in zip(sorted_definitions, row):
                tmp[key_name] = float(val)

            bins.append(tmp)

        sorted_definitions["stat"] = {
            "description": "Uncorrelated statistical uncertainties",
            "treatment": "ADD",
            "type": "UNCORR",
        }
        ret = {"definitions": sorted_definitions, "bins": bins}
        yaml.safe_dump(ret, buffer)

    def export(self, folder_path):
        """Wrapper around export_data and export_uncertainties
        to write both uncertainties and data after filtering to a given folder
        """
        folder_path.mkdir(exist_ok=True)
        # Get the same names as one would use for the filters
        data_path, unc_path = generate_path_filtered_data(folder_path, self.setname)
        # And attach it to the given folder
        data_path = folder_path / data_path.name
        unc_path = folder_path / unc_path.name
        # Export data and uncertainties
        self.export_data(data_path.open("w", encoding="utf-8"))
        self.export_uncertainties(unc_path.open("w", encoding="utf-8"))
        return data_path, unc_path
