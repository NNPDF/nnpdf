"""
This module contains the InconsistentCommonData class which is meant to have all the
methods needed in order to introduce an inconsistency within a Closure Test.
"""

import dataclasses

import pandas as pd

from validphys.coredata import CommonData
from validphys.utils import yaml_safe


@dataclasses.dataclass(eq=False)
class InconsistentCommonData(CommonData):
    """
    Class that inherits all of the methods
    of coredata.CommonData class.

    This class is meant to have all the
    methods needed in order to introduce
    an inconsistency within a Closure Test.
    """

    setname: str
    ndata: int
    commondataproc: str
    nkin: int
    nsys: int
    commondata_table: pd.DataFrame = dataclasses.field(repr=False)
    systype_table: pd.DataFrame = dataclasses.field(repr=False)
    systematics_table: pd.DataFrame = dataclasses.field(default=None, repr=False)
    _systematic_errors: any = dataclasses.field(default=None, init=False)

    @property
    def systematic_errors(self):
        """
        Overrides the systematic_errors method of the CommonData class.

        This is done in order to allow the systematic_errors to be a property
        and hence to be able to assign values to it (setter).
        """
        if self._systematic_errors is None:
            return super().systematic_errors()
        return self._systematic_errors

    @systematic_errors.setter
    def systematic_errors(self, value):
        # Define the setter to allow assignment to systematic_errors
        self._systematic_errors = value

    def select_systype_table_indices(self, treatment_names, names_uncertainties):
        """
        Is used to get the indices of the systype_table that correspond to the
        intersection of the treatment_names and names_uncertainties lists.

        Parameters
        ----------
        treatment_names : list
            list of the names of the treatments that should be selected
            possible values are: MULT, ADD

        names_uncertainties : list
            list of the names of the uncertainties that should be selected
            possible values are: CORR, UNCORR, THEORYCORR, THEORYUNCORR, SPECIAL
            SPECIAL is used for intra-dataset systematics

        Returns
        -------
        systype_tab.index : pd.Index
        """
        # check that names_uncertainties only contains either CORR, UNCORR, THEORYCORR, THEORYUNCORR or SPECIAL
        # if not raise an error
        if not all(
            name in ["CORR", "UNCORR", "THEORYCORR", "THEORYUNCORR", "SPECIAL"]
            for name in names_uncertainties
        ):
            raise ValueError(
                "names_uncertainties should only contain either CORR, UNCORR, THEORYCORR, THEORYUNCORR or SPECIAL"
            )

        # if "SPECIAL", then we need to select the intra-dataset systematics
        if "SPECIAL" in names_uncertainties:
            # avoid circular import error
            from validphys.covmats import INTRA_DATASET_SYS_NAME

            # note: | operator allows to extend the condition so as to also include the names_uncertainties
            systype_tab = self.systype_table[
                (self.systype_table["treatment"].isin(treatment_names))
                & (
                    ~self.systype_table["name"].isin(INTRA_DATASET_SYS_NAME)
                    | self.systype_table["name"].isin(
                        [name for name in names_uncertainties if name != "SPECIAL"]
                    )
                )
            ]

        else:
            systype_tab = self.systype_table[
                (self.systype_table["treatment"].isin(treatment_names))
                & (self.systype_table["name"].isin(names_uncertainties))
            ]

        return systype_tab.index

    def rescale_systematics(self, treatment_names, names_uncertainties, sys_rescaling_factor):
        """
        Rescale the columns of the systematic_errors() that are included in the
        the names_uncertainties list. And return the rescaled table.

        Parameters
        ----------
        treatment_names : list
            list of the names of the treatments that should be rescaled
            possible values are: MULT, ADD

        names_uncertainties : list
            list of the names of the uncertainties that should be rescaled
            possible values are: CORR, UNCORR, THEORYCORR, THEORYUNCORR, SPECIAL
            SPECIAL is used for intra-dataset systematics

        sys_rescaling_factor : float
            factor by which the systematics should be rescaled

        Returns
        -------
        self.systematics_table : pd.DataFrame
        """

        sys_table = self.systematic_errors.copy()

        # select the columns of the systematics_table that should be rescaled
        systype_idx = self.select_systype_table_indices(
            treatment_names=treatment_names, names_uncertainties=names_uncertainties
        )

        # rescale columns of the systematics_table that are included in the index systype_idx
        sys_table.iloc[:, systype_idx - 1] *= sys_rescaling_factor

        return sys_table

    def process_commondata(
        self, treatment_names, names_uncertainties, sys_rescaling_factor, inconsistent_datasets
    ):
        """
        returns a commondata instance
        with modified systematics.
        Note that if commondata.setname
        is not within the inconsistent_datasets or if both ADD and
        MULT are False, then the commondata object
        will not be modified.

        Parameters
        ----------
        treatment_names : list
                            list of the names of the treatments that should be rescaled
                            possible values are: MULT, ADD

        names_uncertainties : list
                            list of the names of the uncertainties that should be rescaled
                            possible values are: CORR, UNCORR, THEORYCORR, THEORYUNCORR, SPECIAL
                            SPECIAL is used for intra-dataset systematics

        sys_rescaling_factor : float, int

        inconsistent_datasets : list
                            list of the datasets for which an inconsistency should be introduced

        Returns
        -------
        validphys.inconsistent_ct.InconsistentCommonData
        """
        new_commondata = self

        if not self.setname in inconsistent_datasets:
            return self

        # needs setter to allow assignment to systematic_errors
        new_commondata.systematic_errors = self.rescale_systematics(
            treatment_names, names_uncertainties, sys_rescaling_factor
        )

        return new_commondata

    def export_uncertainties(self, buffer):
        """
        Same as the export_uncertainties method of the CommonData class.
        The only difference is that systematic_errors is now a property of the class
        and not a method.
        """
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

        for idx, row in self.systematic_errors.iterrows():
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
        yaml_safe.dump(ret, buffer)
