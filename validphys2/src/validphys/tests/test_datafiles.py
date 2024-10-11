"""
    Test all datafiles

    The checks in ``test_all_datasets`` are run for each dataset independently so that one gets one
    failure per dataset in case of problems
"""

import pytest

from validphys.covmats import INTRA_DATASET_SYS_NAME
from validphys.kinematics import xq2map_with_cuts
from validphys.loader import FallbackLoader
from validphys.plotoptions.kintransforms import identity as kintransform_identity

l = FallbackLoader()
all_datasets = sorted(l.implemented_datasets)


def _load_main_and_variants(dataset_name):
    """Given a dataset name, returns a list with the default data and all variants"""
    cds = [l.check_commondata(dataset_name)]
    for variant_name in cds[0].metadata.variants:
        cds.append(l.check_commondata(dataset_name, variant=variant_name))
    return cds


@pytest.mark.parametrize("dataset_name", all_datasets)
def test_all_datasets(dataset_name, data_internal_cuts_new_theory_config):
    """Checks that a dataset can be loaded (together with its variants),
    that the kinematics, uncertainties and data can be read.

    All checks pertaining to a given dataset are done together in this function
    since a broken dataset will likely fail more than one check.
    This avoids polluting the output with many errors for a single dataset.

    Checks:
        1. Loading of data, kinematics, uncertainties
        2. Kinematic coverage is included in the dataframe
        3. A process type is being used (if not legacy)
        4. All variants can be loaded
        5. Uncertainties are either ADD or MULT
        6. The kinematic coverage coverage can be generated
    """
    # Load the data and all its variants
    cds = _load_main_and_variants(dataset_name)
    main_cd = cds[0]

    # Check ndata
    ndata = main_cd.ndata

    # kinematics check
    _ = main_cd.metadata.load_kinematics(drop_minmax=False)
    kin_df = main_cd.metadata.load_kinematics()

    # Check that the kinematic coverage is contained in the kinematics dataframe
    kin_cov = main_cd.metadata.kinematic_coverage
    assert set(kin_cov) <= set(kin_df.columns.get_level_values(0))

    process_type = main_cd.metadata.process_type

    # check whether the kin override is set to the identity
    # and if so, check that the process_type is not simply a string
    kin_override = main_cd.metadata.plotting.kinematics_override
    if isinstance(kin_override, kintransform_identity) and isinstance(process_type, str):
        raise NotImplementedError(f"The {process_type=} is not implemented in process_options")

    elif not isinstance(process_type, str):
        if not process_type.are_accepted_variables(kin_cov):
            raise ValueError(
                f"The dataset {dataset_name} uses {kin_cov} while accepted variables for {process_type} are {process_type.accepted_variables}"
            )

    # load the central data and ndata for every variant
    all_dc = [(cd.metadata.load_data_central(), cd.ndata) for cd in cds]
    # and check they have the same length (it should've been tested internally already)
    assert all(len(i) == ndat for i, ndat in all_dc)

    # check the uncertainties can be loaded
    # note that due to legacy data there might be datasets without data_uncertainties
    # but that would only happen for non-variant (or member 0 of the list)
    # Separate member 0 in that case
    valid_cds = []
    if main_cd.metadata.data_uncertainties:
        valid_cds.append(main_cd)
    valid_cds += cds[1:]

    for cd in valid_cds:
        unc = cd.metadata.load_uncertainties()

        # Check that, if present, the special `stat` key is ADD and UNCORR
        if "stat" in unc:
            stat = unc["stat"]
            assert stat.columns[0][0] == "ADD"
            assert stat.columns[0][1] == "UNCORR"

        intra_dataset = unc.columns.get_level_values("type").isin(
            list(INTRA_DATASET_SYS_NAME) + ["SKIP", "SCALEVAR"]
        )

        # Check that inter dataset correlations are unique
        inter_dataset_corrs = unc.loc[:, ~intra_dataset]
        inter_types = inter_dataset_corrs.columns.get_level_values("type")
        if not inter_types.is_unique:
            raise ValueError(
                f"The inter-dataset uncertainties for {dataset_name} are not unique: {inter_types.value_counts()}"
            )

        # Check that all treatments are either MULT or ADD
        assert set(unc.columns.get_level_values("treatment").unique()) <= {"MULT", "ADD"}

    # Extra checks for non-polarized datasets
    if str(process_type).endswith("_POL"):
        return

    # Legacy datasets with no "new implementation" are skipped
    for cd in valid_cds:
        # check without cuts
        xq2map_with_cuts(cd, False)
