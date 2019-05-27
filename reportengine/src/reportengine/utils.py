"""
Collection of utilities needed in various parts of the library.
"""

import functools
import collections
import pickle
import inspect
import re

#TODO: Support metaclass attributes?
def get_classmembers(cls, *, predicate=None):
    """Return a dictionary mapping member names to their values for a given
    class. Members of base classes will be returned. The dictionary will be
    ordered according to the definition order of the members, starting from
    the most derived classes (in MRO order).

    If predicate is given, only return memvers whose name satisfy it.

    Attributes defined in metaclasses are ignored.
    """
    res = {}
    for base in cls.__mro__:
        for k in base.__dict__:
            if predicate and not predicate(k):
                continue
            if k not in res:
                #See https://bugs.python.org/issue1785
                #Probably not needed without metaclass support.
                try:
                    v = getattr(cls, k)
                except AttributeError:
                    v = base.__dict__[k]
                res[k] = v
    return res

def normalize_name(name):
    """Remove characters not suitable for filenames from the string"""
    bad = re.compile(r'[^\d\w_-]', re.UNICODE)
    return re.sub(bad, '', str(name))


def saturate(func, d):
    """Call a function retrieving the arguments from d by name.
    Variable number of arguments is ignored, and the function cannot
    have positional only arguments."""
    s = inspect.signature(func)
    #This doesn't work as it should
    #ba = s.bind(func, d)
    kwargs = {name: d[name] for name, param in s.parameters.items()
              if param.kind in (param.POSITIONAL_OR_KEYWORD,
                                param.KEYWORD_ONLY)}
    return func(**kwargs)

class comparepartial(functools.partial):
    def __eq__(self, other):
        return (isinstance(other, type(self)) and
                self.func == other.func and self.args == other.args and
                self.keywords == other.keywords
               )

    def __hash__(self):
        return hash(pickle.dumps(self))

class ChainMap(collections.ChainMap):
    def get_where(self, key):
        for i, mapping in enumerate(self.maps):
            try:
                return i, mapping[key]             # can't use 'key in mapping' with defaultdict
            except KeyError:
                pass
        return self.__missing__(key)            # support subclasses that define __missing__

def get_functions(obj):
    """Get the list of members of the object that are functions,
    as an OrderedDict"""
    return collections.OrderedDict(inspect.getmembers(obj, inspect.isfunction))

def get_providers(obj):
    """Return the objects that are likely to make sense as providers.
    This is stricter than what resourceengine uses internally."""
    functions = get_functions(obj)

    def predicate(k,v):
        if hasattr(obj, '__all__'):
            if not k in obj.__all__:
                return False

        return ((not k.startswith('_')) and #hide hidden functions
                inspect.getmodule(v) == obj #require that they are defined in that module
                )

    return collections.OrderedDict((k,v) for k,v in functions.items() if
               predicate(k,v))

def add_highlight(decorator):
    """Add a highlight argument equal to the name of the decorator. This is
    used for helping purposes, to e.g. emphasize which function produces a
    figure."""
    @functools.wraps(decorator)
    def f(*args, **kwargs):
        res = decorator(*args, **kwargs)
        res.highlight = decorator.__name__
        return res
    return f


def ordinal(n):
    """Return an ordinal string for the integer n"""
    residual = n % 10
    teenth = 10 < n % 100 < 20
    if teenth or not residual or residual > 3:
        suffix = 'th'
    else:
        suffix = ('st', 'nd', 'rd')[residual - 1]
    return '%d%s' % (n,suffix)