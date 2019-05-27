# -*- coding: utf-8 -*-
"""
This module provides functionality for manipulationg the
"generalized namespaces" of reportengine.

These work very much like the stack of programming languages like C++ or
Python. The "stack frame" is indexed by a tumple containing the names of the
keys "namespace specification", which can also be a tuple (name, index) for
namespace parts consisting on lists of dictionaries. If two
namespace specifications start with the same sequence, the longer one can see
all the objects of the shorter one.

For example,

.. code-block:: python

    from reportengine import namespaces

    l = [{0: 'x'}, {0:'y'}, {0:'z', 1:'c'},]
    a = {1: 'a'}
    b = {1: 'b'}


    d = {'l': l, 'a': a, 'b': b}
    #d[1] # KeyError
    ns = namespaces.resolve(d, ('a', 'b'))
    ns[1] #b

    ns = namespaces.resolve(d, ('b', 'a'))
    ns[1] #a

    ns = namespaces.resolve(d, ('a' ,'b', ('l', 2)))
    ns[0] #z
    ns[1] #c

Created on Fri Mar  4 15:02:20 2016

@author: Zahari Kassabov
"""
from collections import UserList, UserDict
from collections.abc import Sequence, Mapping

from reportengine.utils import ChainMap, ordinal

__all__ = ('AsNamespace', 'NSList', 'NSItemsDict', 'push_nslevel',
           'expand_fuzzyspec_partial', 'resolve',
           'value_from_spcec_ele')




class AsNamespace:
    def __init__(self, *args, nskey=None, **kwargs):
        self.nskey = nskey
        super().__init__(*args, **kwargs)

    def as_namespace(self):
        return self

    def nsitem(self, item):
        return self[item]

class NSList(AsNamespace, UserList):

    def as_namespace(self):
        return [{self.nskey: item} for item in self]

class NSItemsDict(AsNamespace, UserDict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._nsdicts = {}

    def nsitem(self, item):
        return {self.nskey: self[item]}



class _namespaces: pass

def expand_fuzzyspec_partial(ns, fuzzyspec, currspec=None):
    """Convert a fuzzyspec to a list of specs. Four each key that can't be
    found, yield a tuple:

        key, currspec, ns

    The caller should arange that the key is in the namespace when the
    generator is resumed.
    """
    if not isinstance(ns, ChainMap):
        ns = ChainMap(ns)

    if currspec is None:
        currspec = ()
    if not fuzzyspec:
        return (currspec,)

    ns = resolve(ns, currspec)

    results = []
    #ns = ChainMap(d)
    key, remainder = fuzzyspec[0], fuzzyspec[1:]
    if not key in ns:
        yield key, currspec, ns
    val = ns[key]
    if isinstance(val, Mapping):

        cs_ = (*currspec, key)

        ret = yield from expand_fuzzyspec_partial(ns, remainder, cs_)
        results += [r for r in ret]
    elif isinstance(val, Sequence):
        for i,val_ in enumerate(val):
            if not isinstance(val_, Mapping) and not hasattr(val, 'as_namespace'):
                raise TypeError("Cannot expand non-dict "
                                "list item '%s' (the %s item) of list %s." %
                                (val_, ordinal(i+1), val))
            cs_ = (*currspec, (key, i))

            ret = yield from expand_fuzzyspec_partial(ns, remainder, cs_)
            results += [r for r in ret]
    else:
        raise TypeError("In spec %s, namespace specification '%s' must resolve "
        "to a dict or a list of dicts, not %r." % (currspec,
                                                   key, type(val).__name__))
    return results


def expand_fuzzyspec(ns, fuzzyspec, currspec=None):
    """Return all the nsspecs that spawn from the fuzzyspec.
    Raise ElementNotFound if some part is missing."""
    gen = expand_fuzzyspec_partial(ns, fuzzyspec, currspec)
    try:
        missing, nsspec, _ = gen.send(None)
    except StopIteration as e:
        return e.value
    raise ElementNotFound("Could not resolve a fuzzyspec. "
    "A key is missing: %s, at the level %r." % (missing, nsspec))


def collect_fuzzyspec(ns, key, fuzzyspec, currspec=None):
    """Return the value of key for each spec in the fuzzyspec."""
    specs = expand_fuzzyspec(ns, fuzzyspec, currspec)
    return [resolve(ns, spec)[key] for spec in specs]




def push_nslevel(d, name, value=None):
    if value is None:
        value = {}

    d[name] = value


class ElementNotFound(KeyError): pass

def extract_nsval(ns, ele):

    #Whether the element comes from a shared dictionary that has to be updated
    old = False
    if isinstance(ele, tuple):
        name, index = ele
    else:
        name, index = ele, None

    try:
        if hasattr(ns, 'nsitem'):
            val = ns.nsitem(name)
        else:
            val = ns[name]
            old = True
    except KeyError as e:
        raise ElementNotFound(*e.args   ) from e
    if hasattr(val, 'as_namespace'):
        val = val.as_namespace()
        old = False
    if isinstance(val, Mapping):
        if index is not None:
            raise TypeError("Value %s is a dict, but a "
                "list index was specified" % name)
    elif isinstance(val, Sequence):
        if index is None:
            raise TypeError("Value %s is a list, but no "
            "list index was specified." % name)
        val = val[index]
        if not isinstance(val, Mapping):
            raise TypeError("Value %s in list %s must "
                "be a dictionary, not %s" % (val, ele, type(val)))
    else:
        raise TypeError("Value %s of type %s in %s is not expandable "
                            "as namespace" % (val, type(val), ns))
    if old:
        val = ChainMap({}, val)
    return val



def resolve_partial(ns, spec):
    if not isinstance(ns, ChainMap):
        ns = ChainMap(ns)
    if _namespaces not in ns:
        ns.maps[-1][_namespaces] = {}

    remainder = ()
    if not spec:
        return (), ns
    nsmap = ns[_namespaces]

    if spec in nsmap:
        return tuple(), nsmap[spec]


    for i  in range(len(spec)):

        currspec = spec[:i+1]
        if currspec in nsmap:
            ns = nsmap[currspec]
            continue
        ele = currspec[-1]
        try:
            val = extract_nsval(ns, ele)
        except ElementNotFound:
            #currspec and remainder overlap in one element
            remainder = spec[i:]
            break
        ns = ns.new_child(val)
        nsmap[currspec] = ns

    return remainder, ns


def resolve(d, spec):
    spec = tuple(spec)
    rem, ns = resolve_partial(d, spec)
    if rem:
        raise KeyError("The following parts cannot be expanded %s" % list(rem))
    return ns

def value_from_spcec_ele(ns, ele):
    if isinstance(ele, tuple):
        name, index = ele
        return ns[name][index]
    else:
        return ns[ele]
