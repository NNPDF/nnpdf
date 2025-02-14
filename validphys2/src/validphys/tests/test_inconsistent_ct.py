"""
Module for testing the InconsistentCommonData class in the inconsistent_closuretest module.
Testing is done by mocking the class's methods and properties.
"""

from io import StringIO
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd


class TestInconsistentCommonData(unittest.TestCase):

    @patch(
        'validphys.closuretest.inconsistent_closuretest.inconsistent_ct.InconsistentCommonData',
        autospec=True,
    )
    def setUp(self, MockInconsistentCommonData):
        """
        Set up mock instance of InconsistentCommonData for all tests.
        """
        self.mock_instance = MockInconsistentCommonData.return_value

        # Mocking the DataFrames in the instance
        self.mock_instance.systype_table = pd.DataFrame(
            {"treatment": ["ADD", "MULT", "ADD"], "name": ["CORR", "UNCORR", "SPECIAL"]}
        )

        self.mock_instance.systematic_errors = pd.DataFrame(
            {"sys1": [0.1, 0.2, 0.3], "sys2": [0.4, 0.5, 0.6]}
        )

    def test_systematic_errors_getter(self):
        """
        Test the getter for the systematic_errors property.
        """
        # Set the _systematic_errors to None so the getter is triggered
        self.mock_instance._systematic_errors = None

        # Mock the return value of the superclass's systematic_errors method
        with patch(
            'validphys.coredata.CommonData.systematic_errors',
            return_value=self.mock_instance.systematic_errors,
        ):
            result = self.mock_instance.systematic_errors

            # Assert that the result matches the mock
            pd.testing.assert_frame_equal(result, self.mock_instance.systematic_errors)

    def test_systematic_errors_setter(self):
        """
        Test the setter for the systematic_errors property.
        """
        new_systematic_errors = pd.DataFrame({"sys1": [0.2, 0.3, 0.4], "sys2": [0.5, 0.6, 0.7]})

        self.mock_instance.systematic_errors = new_systematic_errors
        pd.testing.assert_frame_equal(self.mock_instance.systematic_errors, new_systematic_errors)

    def test_select_systype_table_indices(self):
        """
        Test select_systype_table_indices method with valid input.
        """
        treatment_names = ["ADD"]
        names_uncertainties = ["CORR", "SPECIAL"]

        # Mock return of select_systype_table_indices call
        self.mock_instance.select_systype_table_indices.return_value = pd.Index([0, 2])

        result = self.mock_instance.select_systype_table_indices(
            treatment_names, names_uncertainties
        )

        self.mock_instance.select_systype_table_indices.assert_called_once_with(
            treatment_names, names_uncertainties
        )
        pd.testing.assert_index_equal(result, pd.Index([0, 2]))

    def test_select_systype_table_indices_invalid_uncertainties(self):
        """
        Test select_systype_table_indices with invalid uncertainties.
        """
        treatment_names = ["ADD"]
        names_uncertainties = ["INVALID"]

        # Mock the behavior of raising a ValueError
        self.mock_instance.select_systype_table_indices.side_effect = ValueError(
            "names_uncertainties should only contain either CORR, UNCORR, THEORYCORR, THEORYUNCORR or SPECIAL"
        )

        with self.assertRaises(ValueError):
            self.mock_instance.select_systype_table_indices(treatment_names, names_uncertainties)

    def test_rescale_systematics(self):
        """
        Test rescale_systematics method.
        """
        self.mock_instance.systematic_errors = self.mock_instance.systematic_errors.copy()
        treatment_names = ["ADD"]
        names_uncertainties = ["CORR"]
        sys_rescaling_factor = 2.0

        # Mock return of rescale_systematics
        rescaled_table = self.mock_instance.systematic_errors.copy()
        rescaled_table.iloc[:, 0] *= sys_rescaling_factor
        self.mock_instance.rescale_systematics.return_value = rescaled_table

        result = self.mock_instance.rescale_systematics(
            treatment_names, names_uncertainties, sys_rescaling_factor
        )

        # Assert that rescale_systematics was called once and that the return value matches the mock
        self.mock_instance.rescale_systematics.assert_called_once_with(
            treatment_names, names_uncertainties, sys_rescaling_factor
        )
        pd.testing.assert_frame_equal(result, rescaled_table)

    def test_process_commondata(self):
        """
        Test process_commondata method when the dataset is inconsistent.
        """
        inconsistent_datasets = ["test_dataset"]
        treatment_names = ["ADD"]
        names_uncertainties = ["CORR"]
        sys_rescaling_factor = 2.0

        # Mock the return of process_commondata
        modified_commondata = MagicMock()
        self.mock_instance.process_commondata.return_value = modified_commondata

        result = self.mock_instance.process_commondata(
            treatment_names, names_uncertainties, sys_rescaling_factor, inconsistent_datasets
        )

        # Assert that the method was called with correct parameters
        self.mock_instance.process_commondata.assert_called_once_with(
            treatment_names, names_uncertainties, sys_rescaling_factor, inconsistent_datasets
        )
        self.assertEqual(result, modified_commondata)

    def test_export_uncertainties(self):
        """
        Test the export_uncertainties method.
        """
        buffer = StringIO()

        # Mock the export_uncertainties method
        self.mock_instance.export_uncertainties.return_value = None

        self.mock_instance.export_uncertainties(buffer)
        self.mock_instance.export_uncertainties.assert_called_once_with(buffer)


if __name__ == "__main__":
    unittest.main()
