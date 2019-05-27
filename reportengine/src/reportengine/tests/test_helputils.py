import re

from hypothesis import given, assume
from hypothesis.strategies import text

from reportengine import helputils

@given(text(min_size=1, max_size=1000))
def test_sane_wrap(txt):
    assume(txt.strip('\n'))
    ind = '    '
    w = 70
    res = helputils.sane_fill(txt, initial_indent=ind, width=w)
    assert res.startswith(ind)
    assert re.sub(r'\s', '', txt) == re.sub(r'\s', '', res)
    for line in res.splitlines():
        assert line.rfind(' ') < w
