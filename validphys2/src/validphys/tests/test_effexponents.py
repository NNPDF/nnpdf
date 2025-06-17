from validphys.api import API
from validphys.loader import FallbackLoader as Loader
from validphys.scripts.vp_nextfitruncard import PREPROCESSING_LIMS
from validphys.tests.conftest import FIT, FIT_ITERATED
from validphys.utils import yaml_safe


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
    # The runcard of a 2nd iteration fit I ran manually
    # We load it using the context manager because at_input has been modified
    # to load various keys that are not present in the actual runcard for
    # backwards compatibility
    with open(l.check_fit(FIT_ITERATED).path / "filter.yml") as f:
        ite2_runcard = yaml_safe.load(f)

    predicted_ite2_runcard = yaml_safe.load(
        API.iterated_runcard_yaml(fit=FIT, _flmap_np_clip_arg=PREPROCESSING_LIMS)
    )

    # Remove all seed keys from the runcards since these are randomly generated
    # each time, so will not be the same between different runs
    # NB: since some seeds are n3fit-specific, we check that each seed key exists
    runcards = [predicted_ite2_runcard, ite2_runcard]
    fitting_seeds = ["trvlseed", "nnseed", "mcseed"]

    for runcard in runcards:
        for seed in fitting_seeds:
            if seed in runcard:
                runcard.pop(seed)

        # The description can be arbitrarily changed with no effect on the fit
        runcard.pop("description")

    assert predicted_ite2_runcard == ite2_runcard
