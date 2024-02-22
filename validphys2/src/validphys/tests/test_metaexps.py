"""
test_metaexps

Test that the experiments key defined in the commondata meta data, which is
subsequently used for grouping makes sense.

"""
from validphys.api import API
from validphys.loader import Loader


def test_no_systematic_overlaps():
    """Take every available dataset and check that there are no overlapping
    systematics when the grouping is by metadata experiments.

    This is important because we make the assumption that the total covariance
    matrix is block diagonal in metadata experiment, and it is therefore used
    as an optimisation in several places.

    """
    l = Loader()
    ds_names = l.available_datasets
    # TODO: when we have dataset defaults it would be nice to check that none
    # of the default systematics overlap - should be fine for now.
    ds_inputs = [{"dataset": ds} for ds in ds_names]
    res = API.print_systype_overlap(dataset_inputs=ds_inputs, metadata_group="experiment")
    assert isinstance(res, str), f"Overlap found between metadata experiments {res[0]} {res[1]}"
