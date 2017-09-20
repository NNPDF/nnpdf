from hypothesis import given
from hypothesis.strategies import text

from validphys.utils import common_prefix

@given(text(), text())
def test_common_prefix(s1,s2):
    cp = common_prefix(s1,s2)
    assert s1.startswith(cp)
    assert s2.startswith(cp)
    p1 = s1[len(cp):]
    p2 = s2[len(cp):]
    if p1 and p2:
        assert(p1[0] != p2[0])