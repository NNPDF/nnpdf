"""
    Test all datafiles

    The checks in ``test_all_datasets`` are run for each dataset independently so that one gets one
    failure per dataset in case of problems
"""

import pytest

from validphys.covmats import INTRA_DATASET_SYS_NAME
from validphys.kinematics import xq2map_with_cuts
from validphys.loader import FallbackLoader

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
        5. When legacy_data and legacy_theory are included, they "sum up" to legacy
        6. Uncertainties are either ADD or MULT
        7. The kinematic coverage coverage can be generated
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

    # check that the process_type is not simply a string
    if isinstance(process_type, str):
        raise NotImplementedError(f"The {process_type=} is not implemented in process_options")

    # load the central data for every variant
    all_dc = [cd.metadata.load_data_central() for cd in cds]
    # and check they have the same lenght (it should've been tested internally already)
    assert all(len(i) == ndata for i in all_dc)

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

    if "legacy_theory" in main_cd.metadata.variants:
        # First check that if legacy_theory exists by itself, it is because legacy_data also exists
        if "legacy_data" not in main_cd.metadata.variants:
            raise KeyError(
                f"Variant 'legacy_data' must exist whenever 'legacy_theory' does, it doesn't for {dataset_name}"
            )

        # Check that legacy_theory + legacy_data == legacy
        legacy_variant = None
        for cd in valid_cds:
            if cd.metadata.applied_variant == "legacy":
                legacy_variant = cd

        # Now load [legacy_theory, legacy_data]
        legacy_two = l.check_commondata(dataset_name, variant=("legacy_data", "legacy_theory"))

        # The resulting dictionaries of metadata should be equal up to applied_variant
        d1 = dict(legacy_variant.metadata.__dict__)
        d2 = dict(legacy_two.metadata.__dict__)

        d1.pop("applied_variant")
        d2.pop("applied_variant")
        assert (
            d1 == d2
        ), "The tuple of variants (legacy_data, legacy_theory) is not equal to (legacy)"

    # Extra checks for non-polarized datasets
    if str(process_type).endswith("_POL"):
        return

    # Legacy datasets with no "new implementation" are skipped
    for cd in valid_cds:
        # check without cuts
        xq2map_with_cuts(cd, False)
