"""Tests for functions in the validphys.results file.
"""

import numpy as np

from validphys.api import API


def test_groups_central_values_no_table(data_internal_cuts_config):
    """
    Check if the output of `groups_central_values_no_table` agrees with the replica 0 value
    calculated by calling `group_result_table_no_table`. `group_result_table_no_table` also
    computes the predictions for all other replicas.
    """
    compute_only_central = API.groups_central_values_no_table(**data_internal_cuts_config)
    compute_all_replicas = API.group_result_table_no_table(**data_internal_cuts_config)
    np.testing.assert_allclose(compute_only_central, compute_all_replicas['theory_central'])
