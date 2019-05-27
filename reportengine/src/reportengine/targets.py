"""
targets.py

Data structures used to represent computation targets.

Note: At the moment this is in a separate module to resolve circular imports.
Eventually all this should be redesigned away.
"""
from collections import namedtuple

#Target = namedtuple('Target', ('name', 'nsspec', 'extraargs'))
FuzzyTarget = namedtuple('FuzzyTarget', ('name', 'fuzzyspec', 'rootspec', 'extraargs'))