"""
This module contains the InconsistentCommonData class which is meant to have all the
methods needed in order to introduce an inconsistency within a Closure Test.
"""

import dataclasses
from validphys.coredata import CommonData
import pandas as pd


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

    # def systematic_errors(self, central_values=None):
    #     """
    #     Overrides the systematic_errors method of the CommonData class
    #     in order to return the systematics_table attribute.
    #     """
    #     return self.systematics_table

    def select_systype_table_indices(self, treatment_names, names_uncertainties):
        """
        Returns the indices of the systype_table that correspond to the
        names_uncertainties list.

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
            names_uncertainties.remove("SPECIAL")

            # avoid circular import error
            from validphys.covmats import INTRA_DATASET_SYS_NAME

            # note: | operator allows to extend the condition so as to also include the names_uncertainties
            systype_tab = self.systype_table[
                (self.systype_table["treatment"].isin(treatment_names))
                & (
                    ~self.systype_table["name"].isin(INTRA_DATASET_SYS_NAME)
                    | self.systype_table["name"].isin(names_uncertainties)
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
        Rescale the columns of the systematics_table that are included in the
        the names_uncertainties list. And return the rescaled systematics_table

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

        sys_table = self.systematics_table.copy()

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
        new_commondata.systematics_table = self.rescale_systematics(
            treatment_names, names_uncertainties, sys_rescaling_factor
        )

        return new_commondata
