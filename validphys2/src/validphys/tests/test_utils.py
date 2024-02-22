from hypothesis import given
from hypothesis.strategies import lists, text

from validphys.utils import common_prefix


@given(lists(text(), min_size=2))
def test_common_prefix(s):
    cp = common_prefix(*s)
    assert all(i.startswith(cp) for i in s)

    cp = common_prefix(s[0], s[-1])
    p1 = s[0][len(cp) :]
    p2 = s[-1][len(cp) :]
    if p1 and p2:
        assert p1[0] != p2[0]
