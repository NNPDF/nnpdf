"""
checks.py

Checks are preconditions that the arguments of the function must
satisfy at *compile time* in order to execute it. Checks are not
executed by default during the normal operation of the fuction.
Instead they are executed at the time the execution graph is
generated.

Users would typically use the ``make_argcheck`` decorator. This turns
a normal function into a check::


    from reportengine import check, make_argcheck

    @make_argcheck
    def check_argument1(positive_number)
        check(positive_number>1, "The argument should be positive.")

    @check_argument1
    def provider(positive_number:int):
        ...

The checks have almost no effect on the function itself. They only add
append the check to a ``.checks`` attribute of the function object.
However reportengine knows how to look for this attribute and wire the
arguments at compile time.
"""
import functools
from collections.abc import Mapping

from reportengine.utils import saturate
from reportengine.baseexceptions import ErrorWithAlternatives

class CheckError(ErrorWithAlternatives):
    """Error raised by the checking functions"""
    pass

def add_check(f, check):
    """Given a function, ``f`` add a check to it. Users would use
    ``make_argcheck`` or ``make_check`` instead, and this is used to
    implement them."""
    if not hasattr(f, 'checks'):
        f.checks = [check]
    else:
        f.checks.append(check)


def require_one(*args):
    """Ensure that at least one argument is not None."""
    @make_check
    def check(callspec, ns, graph, **kwargs):
        s = set(args)
        in_input_specs = {node.value.resultname for node in graph[callspec].inputs}
        in_ns = {k for k in s if ns.get(k, None) is not None}


        if not (s & (in_ns | in_input_specs)):
            raise CheckError("You need to supply at least one of: %s" % (args,))

    return check

def remove_outer(*args):
    """Set to None all but the innermost values for *args that are not None"""
    @make_check
    def check(ns,**kwargs):
        min_index = len(ns.maps)
        indexes = []
        for arg in args:
            index, val = ns.get_where(arg)
            if val is not None and index < min_index:
                min_index = index
            indexes.append(index)
        for i,arg in zip(indexes, args):
            if i > min_index:
                ns[arg] = None

    return check

def check_positive(var):
    """Ensure that `var` is positive"""
    @make_check
    def check(ns, **kwargs):
        val = ns[var]
        if not val>0:
            raise CheckError(f"'{var}' must be positive, but it is {val!r}.")
    return check


def check_not_empty(var):
    """Ensure that the string ``var`` corresponds to a non empty value in
    the namespace"""
    @make_check
    def check(callspec, ns, graph, **kwargs):
        val = ns[var]
        #Don't just "if val" because we don't know if it's some crazy collection
        if len(val) == 0:
            raise CheckError("'%s' cannot be empty." % var)

    return check

def make_check(check_func):
    """Convert the decorated function in a check. This is mostly
    deprecated, and it is preferred to use ``make_argcheck``. The
    decorated function should take as arguments a mapping, which will
    be a namespace containing the inputs to the provider, and an
    unspecified number of other arguments. The function should raise
    a CheckError in case it fails to validate the arguments (possibly
    using the ``check`` function). The
    effects of mutating the namespace will be visible to the caller.
    """
    @functools.wraps(check_func)
    def decorator(f):
        add_check(f, check_func)
        return f

    return decorator

def make_argcheck(check_func):
    """Convert the decorated function into a check. The inputs should
    have the same name as the inputs of the provider that are to be
    checked. The function should raise a CheckError in case it fails
    to validate the arguments (possibly using the ``check`` function).
    If the function returns a mapping with arguments names askeys, it
    will be used to update the the inputs for the provider. If the
    decorated function doesn't return (retuns ``None``), it will have
    no effect if the checks are succesful."""
    @functools.wraps(check_func)
    @make_check
    def check(ns, *args, **kwargs):
        res = saturate(check_func, ns)
        if res is not None:
            if not isinstance(res, Mapping):
                raise TypeError(f"Bad checking function {check_func}. "
                                "Return value should be None or a mapping, not "
                                f"{type(res)}.")
            ns.update(res)

    return check

def check(cond, *args, **kwargs):
    """Like ``assert`` but not dependent on the interpreter flags, and
    raising a CheckError on failure. ``*args`` and ``**kwargs`` are passed
    to the underlying CheckError."""
    if not cond:
        raise CheckError(*args, **kwargs)
