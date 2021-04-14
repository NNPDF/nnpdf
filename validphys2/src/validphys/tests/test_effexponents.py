import pytest

from validphys.api import API
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT, FIT_ITERATED

from reportengine.compat import yaml


def test_next_runcard():
    """
    Tests that a newly generated iterated fit runcard matches a previously
    generated one (which is stored in FIT_ITERATED).

    Note that since the seeds are generated randomly each time, they will
    be different in the two runcards, so this test does not assert that they
    should be equivalent. In practice this therefore tests that the
    preprocessing exponents and the t0 PDF set are set correctly, and that
    the rest of the information is left untouched.
    """
    l = Loader()
    ite1_fit = l.check_fit(FIT)
    # The runcard of a 2nd iteration fit I ran manually
    ite2_runcard = l.check_fit(FIT_ITERATED).as_input()
    ite2_runcard.pop("pdf")  # Removing the PDF key, it's an artefact of as_input

    # We do this check incase FIT_ITERATED is changed to a new style fit in the
    # future. By checking both namespaces are present, we ensure
    # "dataset_inputs" was added to the fit namespace by the as_input method as
    # opposed to actually being present in the fit runcard.
    if "experiments" in ite2_runcard and 'dataset_inputs' in ite2_runcard:
        # dataset_inputs was added to the as_input for backwards compatibility
        # of the old style fits and wasn't actually present in the fit runcard
        # just like "pdf" above.
        ite2_runcard.pop('dataset_inputs')

    predicted_ite2_runcard = yaml.safe_load(
        API.iterated_runcard_yaml(fit=FIT)
    )

    # Remove all seed keys from the runcards since these are randomly generated
    # each time, so will not be the same between different runs
    # NB: since some seeds are n3fit-specific, we check that each seed key exists
    runcards = [predicted_ite2_runcard, ite2_runcard]
    fitting_seeds = ["seed", "trvlseed", "nnseed", "mcseed"]

    for runcard in runcards:
        if "filterseed" in runcard["closuretest"]:
            runcard["closuretest"].pop("filterseed")
        for seed in fitting_seeds:
            if seed in runcard["fitting"]:
                runcard["fitting"].pop(seed)

    # Check that the actual ite2 runcard matches what vp thinks it should be
    assert predicted_ite2_runcard == ite2_runcard
