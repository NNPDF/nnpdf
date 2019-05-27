# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:31:29 2015

@author: Zahari Kassabov
"""
import inspect
import logging
import functools
import collections
from collections.abc import Mapping
import contextlib
import json


from reportengine.compat import yaml
from reportengine import namespaces
from reportengine.utils import ChainMap, get_classmembers
from reportengine import templateparser
from reportengine.baseexceptions import ErrorWithAlternatives, AsInputError

log = logging.getLogger(__name__)

_config_token = 'parse_'
_produce_token = 'produce_'

def trim_token(s):
    return s.split('_', 1)[1]


class ConfigError(ErrorWithAlternatives):
    pass


class BadInputType(ConfigError, TypeError):
    """Exception that happens when the user enters the wrong input type in the
    config"""
    def __init__(self, param, val, input_type):
        if isinstance(input_type, tuple):
            names = tuple(tp.__name__ for tp in input_type)
        else:
            names = input_type.__name__
        valtype = type(val).__name__
        msg = (f"Bad input type for parameter '{param}': Value '{val}' "
               f"is not of type {names}, but of type '{valtype}'.")
        super().__init__(msg)

class InputNotFoundError(ConfigError, KeyError):
    """Error when the input is not found in the config,"""
    alternatives_header = "Maybe you mistyped %s in one of the following keys?"

def element_of(paramname, elementname=None):
    """Append an elementname and a parentname attribute that will be used
    to generate parsers for lists and mappings of this function."""
    def inner(f):
        nonlocal elementname
        if elementname is None:
            if f.__name__.startswith(_config_token):
                elementname = trim_token(f.__name__)

        f._element_of = paramname
        f._elementname = elementname
        return f
    return inner

def named_element_of(paramname, elementname=None):
    def inner(f):
        element_of(paramname, elementname)(f)
        f._named = True
        return f
    return inner

def _make_element_of(f):
    if getattr(f, '_named', False):
        def parse_func(self, param:dict, **kwargs):
            d = {k: self.trap_or_f(f,  v, f._elementname , **kwargs)
                 for k,v in param.items()}
            return namespaces.NSItemsDict(d, nskey=f._elementname)

        parse_func.__doc__ = "A mapping of %s objects." % f._elementname
    else:
        def parse_func(self, param:list, **kwargs):
            l = [self.trap_or_f(f, elem, f._elementname, **kwargs)
                 for elem in param]
            return namespaces.NSList(l, nskey=f._elementname)
        parse_func.__doc__ = "A list of %s objects." % f._elementname

    #We replicate the same signature for the kwarg parameters, so that we can
    #use that to build the graph.
    list_params = list(inspect.signature(parse_func).parameters.values())[0:2]
    kwarg_params = list(inspect.signature(f).parameters.values())[2:]
    params = [*list_params, *kwarg_params]
    parse_func.__signature__ = inspect.Signature(parameters=params)
    parse_func.__name__ = f'parse_{f._element_of}'
    return parse_func

class ExplicitNode():
    def __init__(self, value):
        self.value = value

def explicit_node(f):
    @functools.wraps(f)
    def f_(*args, **kwargs):
        return ExplicitNode(f(*args,**kwargs))
    return f_


def _parse_func(f):
    """Check that the function has at least one argument, and check that the
    argument corresponds the type declared in the annotation if any."""

    sig = inspect.signature(f)

    try:
        first_param = list(sig.parameters.values())[1]
    except IndexError:
        raise TypeError(("Parser functiom must have at least one "
                        "parameter: %s")
                        % f.__qualname__)

    input_type = first_param.annotation

    @functools.wraps(f)
    def f_(self, val, *args, **kwargs):

        if input_type is not sig.empty:
            if not isinstance(val, input_type):
                raise BadInputType(trim_token(f.__name__), val, input_type)


        return f(self, val, *args, **kwargs)

    return f_

class ElementOfResolver(type):
    """Generate a parsing function for collections of each 'atomic' parsing
    function found in the class, and marked with the relevant decorator."""
    def __new__(cls, name, bases, attrs):
        newattrs = collections.OrderedDict()
        _list_keys = collections.OrderedDict()
        for attr, f in attrs.items():
            if hasattr(f, '_element_of'):
                newattr = _config_token + f._element_of
                if newattr in attrs:
                    raise ValueError("Cannot construct {newattr} from "
                                     "'_element_of' {attr} because it is "
                                     "already declared.")

                #We have to apply parse func in here as well.
                newattrs[newattr] = _make_element_of(_parse_func(f))
                _list_keys[f._element_of] = f._elementname

        newattrs['_list_keys'] = _list_keys

        #We want to make respect the fact that attrs is an ordered dict
        #attrs = {**newattrs, **attrs}
        for k in newattrs:
            attrs[k] = newattrs[k]
        return super().__new__(cls, name, bases, attrs)

class AutoTypeCheck(type):
    """Apply automatically the _parse_func decorator
    to every parsing method fouds in the class."""
    def __new__(cls, name, bases, attrs):
        for k,v in attrs.items():
            if k.startswith(_config_token):
                attrs[k] = _parse_func(v)
        return super().__new__(cls, name, bases, attrs)


class ConfigMetaClass(ElementOfResolver, AutoTypeCheck):
    pass

class Config(metaclass=ConfigMetaClass):

    _traps = ('from_', 'namespaces_')

    def __init__(self, input_params, environment=None):
        if not isinstance(input_params, Mapping):
            raise ConfigError("Failed to process the configuration. Expected "
            "the whole file to resolve to a mapping, but "
            "instead it is %s" % type(input_params))
        self.environment = environment
        self.input_params = input_params

        self._curr_key = None
        self._curr_ns = None
        self._curr_input = None
        self._curr_parents = None
        self._curr_spec = None

        #self.params = self.process_params(input_params)

    @classmethod
    def get_all_parse_functions(cls):
        """Return all defined parse functions, as a dictionary:
        {parsed_element:function}"""
        predicate = lambda x: x.startswith(_config_token)
        return {
            trim_token(k): v
            for (k, v) in get_classmembers(cls, predicate=predicate).items()
        }


    @classmethod
    def get_all_produce_functions(cls):
        """Return all defined parse functions, as a dictionary:
        {parsed_element:function}"""
        predicate = lambda x: x.startswith(_produce_token)
        return {
            trim_token(k): v
            for (k, v) in get_classmembers(cls, predicate=predicate).items()
        }



    def get_parse_func(self, param):
        """Return the function that is defined to parse `param` if it exists.
        Otherwise, return None."""
        func_name = _config_token + param
        try:
            return getattr(self, func_name)
        except AttributeError:
            return None


    def get_produce_func(self, param):
         """Return the function that is defined to produce `param`
         from other inputs.
         Otherwise, return None."""
         func_name = _produce_token + param
         try:
             return getattr(self, func_name)
         except AttributeError:
             return None


    def get_trap_func(self, input_val):
        """If the value has a special meaning that is trapped, return the
        function that handles it. Otherwise, return None"""
        if isinstance(input_val, dict) and len(input_val) == 1:
            k = next(iter(input_val))
            if k in self._traps:
                f = self.get_parse_func(k)
                return functools.partial(f, input_val[k])
        return None


    #TODO: Find out how to avoid duplicaing the functionality
    #of ResourceBuilder.explain_provider
    def explain_param(self, param_name):
        func = self.get_parse_func(param_name)
        start_from = 1

        if func is None:
             func = self.get_produce_func(param_name)
             start_from = 0
             if func is None:
                #TODO: Why not an exception instead of this?
                return None
        result = [func]

        sig = inspect.signature(func)
        for pname, param in list(sig.parameters.items())[start_from:]:
            if self.get_parse_func(pname):
                result.append(('config', pname, self.explain_param(pname)))
            elif self.get_produce_func(pname):
                result.append(('produce', pname, self.explain_param(pname)))
            else:
                result.append(('unknown', pname, param))
        return result



    def trap_or_f(self, f, value, elemname, **kwargs):
        """If the value is a trap, process it, based on the elementname.
        Otherwise just return the result of f(self, value, **kwargs)"""
        tf = self.get_trap_func(value)
        if tf:
            res = tf(elemname, write=False)
            return res[1]
        else:
            return f(self, value, **kwargs)


    def resolve_signature_params(self, f, *, start_from, ns, input_params,
                                 max_index, parents):
        sig = inspect.signature(f)
        kwargs = {}
        put_index = max_index
        for pname, param in list(sig.parameters.items())[start_from:]:
            try:
                index, pval = self.resolve_key(pname,
                                               ns,
                                               input_params= input_params,
                                               max_index=max_index,
                                               parents=parents)
            except KeyError:
                if param.default is not sig.empty:
                    pval = param.default
                    index = max_index
                else:
                    raise

            if index < put_index:
                put_index = index

            kwargs[pname] = pval
        return put_index, kwargs

    @contextlib.contextmanager
    def set_context(self, key=None, ns=None, input_params=None, parents=None,
                    currspec=None):
        """Change the state inside resolve_key.
        Note: This is not for the faint of heart"""
        #Sometimes we just need state, just let's try to not abuse it
        old_key = self._curr_key
        if key is not None:
            self._curr_key = key

        old_ns = self._curr_ns
        if ns is not None:
            self._curr_ns = ns

        old_spec = self._curr_spec
        if currspec is not None:
            self._curr_spec = currspec

        old_inp = self._curr_input
        if input_params is not None:
            self._curr_input = input_params

        old_parents = self._curr_parents
        if parents is not None:
            self._curr_parents = parents

        try:
            yield
        finally:
            self._curr_key = old_key
            self._curr_ns = old_ns
            self._curr_spec = old_spec
            self._curr_input = old_inp
            self._curr_parents = old_parents



    def resolve_key(self, key, ns, input_params=None, parents=None,
                    max_index=None, write=True, currspec=None):
        """Get one key from the input params and put it in the namespace.
        It will be added to the outermost namespace that satisfies all the
        dependencies, but no more levels than `max_index`, if it's given.
        Parents controls the chain of resources that requested this parameter
        (mostly to display errors).
        `write` controls whether the resulting key is to be written to the
        namespace. It only applies to the requested parameter. All dependencies
        will be written to the namespace anyway.
        """
        if parents is None:
            parents = []
        if input_params is None:
            input_params = self.input_params
        with self.set_context(key, ns, input_params, parents, currspec):
            return self._resolve_key(key=key, ns=ns, input_params=input_params,
                parents=parents, max_index=max_index, write=write)


    def _value_and_line_info(self, key, input_params):
        if isinstance(input_params, ChainMap):
            ind, val = input_params.get_where(key)
            inp = input_params.maps[ind]
        else:
            inp = input_params
            val = input_params[key]
        if hasattr(inp, 'lc'):
            lineinfo = inp.lc.item(key)
        else:
            lineinfo = None
        return val, lineinfo


    def _resolve_key(self, key, ns, input_params=None, parents=None,
                    max_index=None, write=True):
        """Do not call this directly. Use resolve_key instead"""



        #TODO: Is there any way to break this down into pieces?
        #      Maybe move handling of produce rules to resourcebuilder?

        if max_index is None:
            max_index = len(ns.maps) -1


        nsindex = nsval = None
        finindex = finval = None
        if key in ns:
            ind, val = ns.get_where(key)
            if ind <= max_index:
                nsindex, nsval = ind, val
            finindex, finval = ind, val


        newparents=[*parents, key]
        produce_func = self.get_produce_func(key)

        #TODO: This logic is pretty horrible...

        if not key in input_params:
            if produce_func:
                try:
                    put_index, kwargs = self.resolve_signature_params(produce_func,
                                            start_from=0 ,ns=ns,
                                            parents=newparents,
                                            input_params=input_params,
                                            max_index=max_index)
                except InputNotFoundError as e:
                    log.debug(f"Can't satisfy production rule for {key}")
                    if nsindex is not None:
                        return nsindex, nsval
                    raise e
                if nsindex is not None and nsindex <= put_index:
                    return nsindex, nsval
                val = produce_func(**kwargs)
                if write:
                    self._write_val(ns, key, val, put_index)
                return put_index, val

            if finindex is not None:
                return finindex, finval
            msg = "A parameter is required: {key}.".format(key=key)
            if parents:
                msg += "\nThis is needed to process:\n"
                msg += '\ntrough:\n'.join(' - ' + str(p) for
                                          p in reversed(parents))
            #alternatives_text = "Note: The following similarly spelled "
            #                     "params exist in the input:"

            raise InputNotFoundError(msg, key, alternatives=input_params.keys())

        put_index = max_index
        input_val, lineinfo = self._value_and_line_info(key, input_params)

        trap_func = self.get_trap_func(input_val)

        if trap_func:
            #TODO: Think about this interface
            return trap_func(key)


        f = self.get_parse_func(key)
        if f:
            put_index, kwargs = self.resolve_signature_params(f,
                                    start_from=1 ,ns=ns,
                                    parents=newparents,
                                    input_params=input_params,
                                    max_index=max_index)

            if nsindex is not None and nsindex <= put_index:
                return nsindex, nsval
            try:
                val = f(input_val, **kwargs)
            except ConfigError as e:
                if lineinfo:
                    log.error(f"Failed processing key '{key}' at line "
                              f"{lineinfo[0]+1}, position {lineinfo[1]+1}.")
                else:
                    log.error(f"Failed processing key {key}.")
                raise e
        elif nsindex is not None:
            return nsindex, nsval
        elif produce_func:
            raise ConfigError(("The key '%s' confilcts with a production rule "
            "and no parser is defined.") % key)
        else:
            self.check_key_is_allowed(key, input_val)


            #Recursively parse dicts
            if isinstance(input_val, dict):
                put_index = 0
                val = {}
                res_ns = ns.new_child(val)
                inputs = ChainMap(input_val, input_params)
                for k in input_val.keys():
                    self.resolve_key(k, res_ns, inputs,
                                     parents=[*parents, key],
                                     max_index = 0
                                    )
            #Recursively parse lists of dicts
            elif (isinstance(input_val, list) and
                 all(isinstance(x, dict) for x in input_val)):
                put_index = 0
                val = []
                for linp in input_val:
                    lval = {}
                    res_ns = ns.new_child(lval)
                    inputs = ChainMap(linp, input_params)
                    for k in linp.keys():
                        self.resolve_key(k, res_ns, inputs,
                                         parents=[*parents, key],
                                         max_index = 0
                                        )
                    val.append(lval)


            else:
                val = input_val
        if write:
            self._write_val(ns, key, val, put_index)


        return put_index, val

    def _write_val(self, ns, key ,val, put_index):
        #TODO: Need to fix this better
        if isinstance(val, ExplicitNode):
            ns[key]=val
        else:
            ns.maps[put_index][key] = val

    def process_fuzzyspec(self, fuzzy, ns, parents=None, initial_spec=None):
        if parents is None:
            parents = []
        gen = namespaces.expand_fuzzyspec_partial(ns, fuzzy, currspec=initial_spec)
        while True:
            try:
                key, currspec, currns = next(gen)
            except StopIteration as e:
                return e.value
            except TypeError as e:
                raise ConfigError("Error when processing namespace "
                "specification %s: %s" % (fuzzy, e))
            else:
                self.resolve_key(key, currns, parents=[*parents, currspec],
                                 currspec=currspec)

    def process_all_params(self, input_params=None, *,ns=None):
        """Simple shortcut to process all paams in a simple namespace, if
        possible."""
        if input_params is None:
            input_params = self.input_params

        if ns is None:
            ns = ChainMap()
        for param in input_params:
            if param not in ns:
                self.resolve_key(param, ns, input_params=input_params)
        return ns

    def _parse_actions_gen(self, actions, currspec=()):
        if isinstance(actions, dict):
            for k,v in actions.items():
                yield from self._parse_actions_gen(v, (*currspec, k))
        elif isinstance(actions, list):
            for v in actions:
                if isinstance(v, dict):
                    if len(v) != 1:
                        raise ConfigError(("Invalid action specification %s. "
                        "Must be a scalar or a mapping with exactly one key") % v)
                    k = next(iter(v.keys()))
                    args = v[k]
                    if not isinstance(args, dict):
                        raise ConfigError("Action arguments must be "
                        "a mapping if present: %s" % k)
                    yield k, currspec, () ,tuple(args.items())
                elif isinstance(v, str):
                    yield v, currspec, (), ()
                else:
                    raise ConfigError("Unrecognized format for actions. "
                    "Must be a string or mapping, not '%s'" %v)
        else:
            raise ConfigError("Unrecognized format for actions")

    def check_key_is_allowed(self, key, value):
        if not hasattr(self, 'allowed_keys'):
            return
        if key not in self.allowed_keys:
            msg = "Keyword '%s' is not allowed in the configuration."%key
            alternatives = self.allowed_keys.keys() | self.get_all_parse_functions().keys()
            raise ConfigError(msg, key, alternatives)
        good_tp = self.allowed_keys[key]
        if not isinstance(value, good_tp):
            raise BadInputType(key, value, good_tp)


    def _rewrite_old_actions(self, res):
        rewrite_pieces = []
        for act in res:
            name, fuzzyspec, root, extraargs = act
            assert not root
            if extraargs:
                eapieces = []
                for ea in extraargs:
                    k, v = ea
                    eapiece = f'{k}={json.dumps(v)}'
                    eapieces.append(eapiece)
                eapart = f"({', '.join(eapieces)})"
            else:
                eapart = ''

            rewrite_pieces.append(f"  - {'::'.join(fuzzyspec)}{' ' if fuzzyspec else ''}{name}{eapart}")
        allactions = '\n'.join(rewrite_pieces)
        return f"\nactions_:\n{allactions}\n"



    def parse_actions_(self, actions:list):
        """A list of action specifications. See documentation and examples for details."""
        if all(isinstance(act,str) for act in actions):
            try:
                return [templateparser.string_to_target(act) for act in actions]
            except templateparser.BadTemplate as e:
                raise ConfigError(e) from e
        allacts = [list(self._parse_actions_gen(act)) for act in actions]
        res = [act for acts in allacts for act in acts]
        newactions = self._rewrite_old_actions(res)
        log.warning(f"Old style actions_ are deprecated. Rewrite them like this:\n{newactions}")
        return res

    def parse_namespaces_(self, value:str, element):
        """Convert the referenced key into a list of namespaces.
        The value is interpreted in the same way as the content of the
        {@with@} blocks in the report template. That is, the namespace
        elements are delimited by ``::``."""
        if not isinstance(value, str):
            raise ConfigError(f'Argument of namespace_ must '
                              'be str, not {type(element)}')

        currspec = self._curr_spec
        ns = self._curr_ns.maps[-1]


        fuzzy = templateparser.tokenize_fuzzy(value)
        log.debug(f"Obtaining namespaces from specs {fuzzy}.")


        if currspec:
            raise ConfigError("namespaces_ is only supported as a top "
                f"level namespace, specificaation, but {fuzzy} it is being "
                f"resolved inside {currspec}.")


        specs = self.process_fuzzyspec(fuzzy, ns, self._curr_parents, initial_spec=currspec)

        res = [namespaces.resolve(ns, spec) for spec in specs]
        ns[element] = res

        return 0, res


    #TODO: This interface is absolutely horrible, but we need to do a few
    #more examples (like 'zip_') to see how it generalizes.
    def parse_from_(self, value:(str, type(None)), element, write=True):
        """Take the key from the referenced element,
        which should be another input resource  and resolve as a dict.
        If the value is None, retrieve the key from the current namespace.
        """

        ns = self._curr_ns
        input_params = self._curr_input
        parents = self._curr_parents

        if value is None:
            return self.resolve_key(element, ns, input_params=input_params,
                                    parents=parents, write=write)



        nokey_message = ("Could retrieve element %s from namespace. "
                          "No such key" %
                              (element,))

        #The write=False is needed here because of the general stupidy of the
        #framework. Because there is no way to tell how deep in the inputs
        #we found the key, it may end up completely misplaced. An upgrade
        #of the code should fix this.
        max_index, tip = self.resolve_key(value, ns, input_params=input_params,
                                        parents=parents, write=False)

        if hasattr(tip, 'as_input'):
            try:
                d = tip.as_input()
            except AsInputError as e:
                raise ConfigError("Could not process '%s' as input: %s" % (value, e)) from e
            try:
                ele_input = d[element]
            except KeyError as e:
                raise ConfigError(nokey_message) from e

            new_input = ChainMap({element:ele_input},input_params)
            #This is done due to the (crazy) put_index logic
            new_ns = ns.new_child()
            new_ns[value] = tip

            ind,val = self.resolve_key(element, new_ns,
                input_params=new_input,
                parents=parents, write=write, max_index=max_index)
            ns.update(new_ns)
            return ind,val


        elif isinstance(tip, dict):
            d = tip
            try:
                res  = d[element]
            except KeyError as e:
                raise ConfigError(nokey_message) from e
            if write:
                ns[element] = res
            return 0, res
        else:
            bad_message = ("Unrecognized value for from_: "
            "The value must resolve to a mapping or an object implementing "
            "as_input")
            raise ConfigError(bad_message)

    def __getitem__(self, item):
        return self.input_params[item]

    def __iter__(self):
        return iter(self.input_params)

    def __len__(self):
        return len(self.input_params)

    def __contains__(self, item):
        return item in self.input_params

    @classmethod
    def from_yaml(cls, o, *args, **kwargs):
        try:
            return cls(yaml.round_trip_load(o), *args, **kwargs)
        except yaml.error.YAMLError as e:
            raise ConfigError(f"Failed to parse yaml file: {e}")
