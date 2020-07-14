import pytest

from validphys.api import API
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT, FIT_ITERATED

from reportengine.compat import yaml


def test_next_runcard():
    l = Loader()
    ite1_fit = l.check_fit(FIT)
    # The runcard of a 2nd iteration fit I ran manually
    ite2_runcard = l.check_fit(FIT_ITERATED).as_input()
    ite2_runcard.pop("pdf")  # Removing the PDF key, it's an artefact of as_input

    predicted_ite2_runcard = yaml.safe_load(
        API.next_effective_exponents_yaml(fit=FIT)
    )
    # Check that the actual ite2 runcard matches what vp thinks it should be
    assert predicted_ite2_runcard == ite2_runcard
