"""
test_theorycovariance.py

Tests for the theory covariance matrices built in
:py:mod:`validphys.theorycovariance.construction`, and in particular for the
``nnfit_theory_covmat`` production rule, which is what ``vp-setupfit`` asks for
whenever a fit runcard contains a ``theorycovmatconfig``.
"""

import numpy as np
import pytest

from reportengine.table import savetable
from validphys.api import API

# The user covmat is looked up with ``Loader.check_vp_output_file``, which
# rejects absolute paths, so it has to be found relative to the cwd. The name is
# made distinctive so that it cannot collide with anything in the vp cache.
USER_COVMAT_FILENAME = "test_user_covmat.csv"
POINT_PRESCRIPTION = "3 point"


@pytest.fixture(scope="module")
def user_covmat_on_disk(tmp_path_factory, thcovmat_config):
    """Write out a "user" theory covmat and return ``(directory, matrix)``.

    The matrix written out is the ``3 point`` scale variation covmat for the very
    same data, so that once it has been read back by ``fromfile_covmat`` the
    resulting ``user_covmat`` is numerically identical to ``theory_covmat_custom``.
    This exercises the real cut and index handling of ``fromfile_covmat`` and
    gives the tests below an expected value with no hardcoded numbers in it.

    ``theory_covmat_custom`` is float32; it is promoted to float64 before being
    written, because the shortest repr of a float32 only round trips to float32
    and would introduce a ~1e-8 relative error.
    """
    covmat = API.theory_covmat_custom(
        point_prescriptions=[POINT_PRESCRIPTION], **thcovmat_config
    ).astype(np.float64)
    path = tmp_path_factory.mktemp("user_covmat")
    savetable(covmat, path / USER_COVMAT_FILENAME)
    return path, covmat


@pytest.mark.parametrize(
    ("covmat_config", "expected_factor"),
    [
        # Scale variations *and* a user covmat. This is the branch that used to
        # fail during graph resolution with
        # ``KeyError: 'group_dataset_inputs_by_metadata'``.
        (
            {"point_prescriptions": [POINT_PRESCRIPTION], "user_covmat_path": USER_COVMAT_FILENAME},
            2.0,
        ),
        # Scale variations only.
        ({"point_prescriptions": [POINT_PRESCRIPTION]}, 1.0),
        # User covmat only.
        ({"user_covmat_path": USER_COVMAT_FILENAME}, 1.0),
        # User covmat only, rescaled by ``mult_factor``.
        ({"user_covmat_path": USER_COVMAT_FILENAME, "mult_factor": 2.5}, 2.5),
    ],
    ids=["scalevar_and_user", "scalevar_only", "user_only", "user_only_rescaled"],
)
def test_nnfit_theory_covmat(
    monkeypatch, thcovmat_config, user_covmat_on_disk, covmat_config, expected_factor
):
    """Every branch of ``produce_nnfit_theory_covmat`` must resolve and return a
    covmat indexed in runcard order on both axes.

    The ``scalevar_and_user`` case is a regression test: combining
    ``point_prescriptions`` with ``user_covmat_path`` used to raise
    ``KeyError: 'group_dataset_inputs_by_metadata'`` while the reportengine graph
    was being built, because ``total_theory_covmat`` carried a check whose
    arguments are not resolvable in the ``theorycovmatconfig`` namespace.
    """
    covmat_dir, scalevar_covmat = user_covmat_on_disk
    monkeypatch.chdir(covmat_dir)

    # Resolving and executing this at all is the actual regression.
    covmat = API.nnfit_theory_covmat(**covmat_config, **thcovmat_config)

    runcard_index = API.data_index(**thcovmat_config).droplevel(0)
    process_index = API.procs_index(**thcovmat_config).droplevel(0)

    # Guard against this test quietly becoming vacuous. It only has teeth as long
    # as grouping DATA_THCOVMAT by process actually reorders it: LHCB_Z0_8TEV_MUON_Y
    # is DY NC and so joins the group opened by CMS_Z0J_8TEV_PT-Y, overtaking
    # ATLAS_WJ_8TEV_WP-PT, which is DY CC. If DATA_THCOVMAT is ever changed into an
    # order that the grouping leaves alone, this test must be given its own
    # dataset_inputs rather than being left silently trivial.
    assert not runcard_index.equals(process_index), (
        "Grouping DATA_THCOVMAT by process no longer reorders it, so the "
        "assertions below would hold even if the covmat were written out in "
        "process order instead of runcard order."
    )

    assert covmat.index.droplevel(0).equals(runcard_index)
    assert covmat.columns.equals(covmat.index)

    expected = scalevar_covmat.reindex(index=covmat.index, columns=covmat.columns)
    np.testing.assert_allclose(covmat.to_numpy(), expected_factor * expected.to_numpy(), rtol=1e-8)
