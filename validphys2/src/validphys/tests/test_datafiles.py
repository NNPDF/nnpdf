"""
    Test all datafiles

    All tests are under ``test_all_datasets`` ran with all datasets so that one gets one failure per dataset in case of problems
"""

import pytest

from validphys.covmats import INTRA_DATASET_SYS_NAME
from validphys.loader import FallbackLoader
from validphys.plotoptions.kintransforms import identity as kintransform_identity

l = FallbackLoader()
all_datasets = sorted(l.implemented_datasets)


def _load_main_and_variants(dataset_name):
    """Given a dataset name, returns a list with the default load and all variants"""
    cds = [l.check_commondata(dataset_name)]
    for variant_name in cds[0].metadata.variants:
        cds.append(l.check_commondata(dataset_name, variant=variant_name))
    return cds


@pytest.mark.parametrize("dataset_name", all_datasets)
def test_all_datasets(dataset_name):
    """checks that a dataset can be loaded (together with its variants),
    that the kinematics, uncertainties and data can be read
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
        # Skip for the time being the processes for which there is no implementation but have been
        # merged to master: issue #1991
        if process_type not in ("HQP_MQQ", "INC"):
            raise NotImplementedError(f"The {process_type=} is not implemented in process_options")

    elif not isinstance(process_type, str):
        if not process_type.are_accepted_variables(kin_cov):
            raise ValueError(
                f"The dataset {dataset_name} uses {kin_cov} while accepted variables for {process_type} are {process_type.accepted_variables}"
            )

    # load the central data for every variant
    all_dc = [cd.metadata.load_data_central() for cd in cds]
    # and check they have the same lenght (it should've been tested internally already)
    assert all(len(i) == ndata for i in all_dc)

    # check the uncertainties can be loaded
    # note that due to legacy data there might be datasets without data_uncertainties
    # but that would only happen for non-variant (or member 0 of the list)
    all_unc = [cd.metadata.load_uncertainties() for cd in cds[1:]]
    if main_cd.metadata.data_uncertainties:
        all_unc.insert(0, main_cd.metadata.load_uncertainties())

    for unc in all_unc:
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
