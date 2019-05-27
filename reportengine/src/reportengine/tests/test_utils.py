"""
Tests for the utils module.
"""

from reportengine.utils import get_classmembers

class Meta(type):
    @property
    def foo(self):
        return "foo"

class A(metaclass=Meta):
    b = 1
    c = 2

class B(A):
    b=4
    d=1

def test_get_classmembers():
    m = get_classmembers(B, predicate=lambda x: not x.startswith('_'))
    assert list(m) == ['b', 'd', 'c']
    assert m['b'] == 4

