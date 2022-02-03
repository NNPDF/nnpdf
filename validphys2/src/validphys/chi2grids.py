"""
chi2grids.py

Compute and store χ² data from replicas, possibly keeping the correlations
between pseudorreplica fluctuations between different fits. This is applied
here to parameter determinations such as those of αs.
"""
import logging
from collections import namedtuple

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table
from validphys.calcutils import calc_chi2

PseudoReplicaExpChi2Data = namedtuple(
    "PseudoReplicaChi2Data", ["group", "ndata", "chi2", "nnfit_index"]
)


log = logging.getLogger(__name__)


def computed_pseudoreplicas_chi2(
    fitted_make_replicas,
    fitted_replica_indexes,
    group_result_table_no_table,  # to get the results already in the form of a dataframe
    groups_sqrtcovmat,
):
    """Return a dataframe with the chi² of each replica with its corresponding
    pseudodata (i.e. the one it was fitted with). The chi² is computed by group.
    The index of the output dataframe is
        ``['group',  'ndata' , 'nnfit_index']``
    where ``nnftix_index`` is the name of the corresponding replica
    """
    # Stack the replica pseudodata to have the prediction shape
    r_data = np.stack(fitted_make_replicas[0], axis=1)

    # Drop data central and theory central which is not useful here
    r_prediction = group_result_table_no_table.drop(columns=["data_central", "theory_central"])

    # Now compute the chi2 in a per-group basis
    diff = r_prediction - r_data
    group_level = r_prediction.index.get_level_values("group")

    # Save the results in a dataframe similar (but not equal) to the old one
    df_output = []
    for group in group_level.unique():
        group_diff = diff.loc[group_level == group]
        its_covmat = groups_sqrtcovmat[group_level == group][group]
        chi2_per_replica = calc_chi2(its_covmat, group_diff)
        ndata = len(group_diff)
        for i, chi2 in zip(fitted_replica_indexes, chi2_per_replica):
            df_output.append(PseudoReplicaExpChi2Data(group, ndata, chi2, i))

    df = pd.DataFrame(df_output, columns=PseudoReplicaExpChi2Data._fields)
    df.set_index(["group", "ndata", "nnfit_index"], inplace=True)
    df.sort_index(inplace=True)
    return df


fits_computed_pseudoreplicas_chi2 = collect(
    computed_pseudoreplicas_chi2, ("fits", "fitenvironment", "fitpdf")
)

dataspecs_computed_pseudorreplicas_chi2 = collect(computed_pseudoreplicas_chi2, ("dataspecs",))

@table
def export_computed_pseudoreplicas_chi2(computed_pseudoreplicas_chi2):
    return computed_pseudoreplicas_chi2


@table
def export_fits_computed_pseudoreplicas_chi2(fits_computed_pseudoreplicas_chi2):
    """Hack to force writting the CSV output"""
    return fits_computed_pseudoreplicas_chi2
