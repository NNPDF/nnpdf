# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 21:18:06 2015

@author: zah
"""
from __future__ import generator_stop

from collections import namedtuple, defaultdict, OrderedDict
from collections.abc import Sequence
import os
import logging
import inspect
import functools
from abc import ABCMeta
import warnings

import curio

from reportengine import dag
from reportengine import namespaces
from reportengine.configparser import InputNotFoundError, BadInputType, ExplicitNode
from reportengine.checks import CheckError
from reportengine.utils import ChainMap
from reportengine.targets import FuzzyTarget

log = logging.getLogger(__name__)

RESOURCE = "resource"
PROVIDER = "provider"

#These represent the final actions we are interested in executiong.



EMPTY = inspect.Signature.empty

class resultkey:pass

#TODO: Fix help system for collect
class collect:
    resultkey = resultkey
    def __init__(self, function, fuzzyspec, *,element_default=EMPTY):
        self.fuzzyspec = fuzzyspec
        if isinstance(function, collect):
            raise TypeError("Unsupported collect call taking a colect instance as first argument. "
                            "Pass a string with the name of the variable instead.")
        if not isinstance(function, str) and not hasattr(function, '__name__'):
            raise TypeError("Invalid argument 'function'. Must be either a "
                        f"named function or a string, not {type(function)}")
        self.function = function
        self.element_default = element_default

    def __call__(self, ns, nsspec):
        ns = namespaces.resolve(ns, nsspec)
        return list(ns[resultkey].values())

class target_map(collect):
    class targetlenskey:pass
    def __init__(self, targets):
        self._targets = targets

    def get_targets(self, nsspec):
        return [target._replace(rootspec=nsspec) for target in self._targets]


    def __call__(self, ns, nsspec):
        ns = namespaces.resolve(ns,nsspec)
        return ns[resultkey]

def add_to_dict_flag(result, ns ,origin, target, index):
    log.debug("Setting element %s of %r from %r", index, target, origin)
    namespaces.resolve(ns, target.nsspec)[collect.resultkey][index] = result

async def _async_identity(f, *args, **kwargs):
    return f(*args, **kwargs)

class provider:
    """Decorator intended to be used for the functions that are to
    be exposed as providers, either directly or trough more specialized
    decorators in reportengine."""
    def __init__(self, f):
        functools.update_wrapper(self, f)
        self.f = f

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)

class Node(metaclass = ABCMeta):
    pass

CallSpec = namedtuple('CallSpec', ('function', 'kwargs', 'resultname',
                                   'nsspec'))

CollectSpec = namedtuple('CollectSpec', ('function', 'kwargs', 'resultname',
                                   'nsspec'))

CollectMapSpec = namedtuple('CollectMapSpec', ('function', 'kwargs' ,'resultname', 'nsspec'))

Node.register(CallSpec)
Node.register(CollectSpec)
Node.register(CollectMapSpec)


#TODO; Improve namespace spec
def print_callspec(spec, nsname = None):

    if nsname is None:
        res = spec.resultname
    else:
        res = "nsname[{!r}]".format(spec.resultname)
    callargs = ', '.join("%s=%s"% (kw, kw) for kw in spec.kwargs)
    try:
        f = spec.function.__qualname__

    #functools.partial annoyingly doesn't wrap
    except AttributeError:
        f = spec.function.func.__qualname__
        callargs += ", " + ", ".join("%s=%s" % (kw,val) for
                                     kw,val in spec.function.keywords.items())

    return "{res} = {f}({callargs})".format(res=res,
                                        f=f,
                                        callargs=callargs)

CallSpec.__str__ = print_callspec

def check_types(f, ns):
    s = inspect.signature(f)
    for param_name, param_value in s.parameters.items():
        tps = param_value.annotation
        if tps is not s.empty and param_name in ns:
            val = ns[param_name]
            if not isinstance(val, tps):
                raise BadInputType(param_name, val, tps)

class ResourceExecutor():

    def __init__(self, graph, rootns, environment=None, perform_final=True):
        self.graph = graph
        self.rootns = rootns
        self.environment = environment
        self._node_flags = defaultdict(lambda: set())
        self.perform_final = perform_final

    def resolve_callargs(self, callspec):
        function, kwargs, resultname, nsspec = callspec
        namespace = namespaces.resolve(self.rootns, nsspec)
        kwdict = {kw: namespace[kw] for kw in kwargs}
        if hasattr(function, 'prepare') and self.perform_final:
            prepare_args = function.prepare(spec=callspec,
                                            namespace=self.rootns,
                                            environment=self.environment,)
        else:
            prepare_args = {}

        return kwdict, prepare_args

    def execute_sequential(self):
        for node in self.graph:
            callspec = node.value
            if isinstance(callspec, (CollectSpec, CollectMapSpec)):
                #my_ns = namespaces.resolve(self.rootns, callspec.nsspec)
                result = callspec.function(self.rootns, callspec.nsspec)
            else:
                result = self.get_result(callspec.function,
                                         *self.resolve_callargs(callspec),
                                         perform_final=self.perform_final)
            self.set_result(result, callspec)

    #This needs to be a staticmethod, because otherwise we have to serialize
    #the whole self object when passing to multiprocessing.
    @staticmethod
    def get_result(function, kwdict, prepare_args, perform_final=True):
        fres =  function(**kwdict)
        if hasattr(function, 'final_action') and perform_final:
            return function.final_action(fres, **prepare_args)
        return fres

    def set_result(self, result, spec):
        function, _, resultname, nsspec = spec
        namespace = namespaces.resolve(self.rootns, nsspec)
        put_map = namespace.maps[1]
        log.debug("Setting result for %s %s", spec, nsspec)

        put_map[resultname] = result

        for action, args in self._node_flags[spec]:
            action(result, self.rootns ,spec, **dict(args))



    async def _run_parallel(self, deps, completed_spec):
        try:
            runnable_specs = deps.send(completed_spec)
        except StopIteration:
            return
        pending_tasks = {}
        tg = curio.TaskGroup()
        for pending_spec in runnable_specs:
            if isinstance(pending_spec, (CollectSpec, CollectMapSpec)):
                remote_coro = _async_identity(pending_spec.function,
                          self.rootns, pending_spec.nsspec)
            else:
                remote_coro = curio.run_in_process(self.get_result,
                                       pending_spec.function,
                                       *self.resolve_callargs(pending_spec))
            pending_task = await tg.spawn(remote_coro)
            pending_tasks[pending_task] = pending_spec

        next_runs =  curio.TaskGroup()

        async for completed_task in tg:
            try:
                result = await completed_task.join()
            except curio.TaskError as e:
                raise curio.KernelExit() from e

            new_completed_spec = pending_tasks.pop(completed_task)
            self.set_result(result, new_completed_spec)
            next_runs_coro = self._run_parallel(deps, new_completed_spec)
            await next_runs.spawn(next_runs_coro)
        async for t in next_runs:
            await t.join()





    def execute_parallel(self):
        if 'MAX_WORKER_PROCESSES' in os.environ:
            try:
                nworkers = int(os.environ['MAX_WORKER_PROCESSES'])
            except Exception as e:
                log.warning('Could not interpret the value of the environment variable MAX_WORKER_PROCESSES as an integer. Ignoring it.')
            if nworkers>=1:
                log.debug(f"Setting MAX_WORKER_PROCESSES to {nworkers} from the environment.")
                curio.workers.MAX_WORKER_PROCESSES = nworkers
            else:
                log.warning('The environment variable MAX_WORKER_PROCESSES must be greater >=1. Ignoring it.')

        deps = self.graph.dependency_resolver()
        #https://github.com/dabeaz/curio/issues/72
        kernel = curio.Kernel()
        async def main():
            await self._run_parallel(deps, None)

        with kernel:
            try:
                kernel.run(main())
            except KeyboardInterrupt:
                log.info("Canceling remaining tasks")
                raise
            except curio.KernelExit as e:
                raise e.__cause__

    def __str__(self):
        return "\n".join(print_callspec(node.value) for node in self.graph)


class ResourceError(Exception):
    def __init__(self, name, message, parents):
        self.name = name
        self.message = message
        if not parents:
            parents = ('Target specification',)
        self.parents = parents

    def __str__(self):
        return "Could not process the resource '%s', required by:\n%s\n%s"%(
                self.name, '\n'.join(' - ' + p for p in self.parents),
                self.message)

class ResourceNotUnderstood(ResourceError, TypeError): pass


class ResourceBuilder(ResourceExecutor):

    def __init__(self, input_parser, providers, fuzzytargets, environment=None, perform_final=True):
        self.input_parser = input_parser

        if not isinstance(providers, Sequence):
            providers = [providers]
        self.providers = providers
        self.fuzzytargets = fuzzytargets

        rootns = ChainMap()
        graph = dag.DAG()
        super().__init__(graph, rootns, environment, perform_final)


    def is_provider_func(self, name):
        return any(hasattr(provider, name) for provider in self.providers)

    def get_provider_func(self, name):
        for provider in self.providers:
            func = getattr(provider, name, False)
            if func:
                return func
        raise AttributeError("No such provider function: %s" % name)

    def explain_provider(self, provider_name):
        """Collect the dependencies of a given provider name, by analizing
        the function signature.
        The result is a list where the first element is the provider function.
        The next elements are tuples of the form ``(type, name, value)`` where

         - ``type`` is one of ``('provider', 'config', 'produce', 'unknown')``
         - ``name`` is the name of the dpenendncy
         - ``value`` depends on the type:
          - For providers, it is the output of ``self.explain_provider(name)``.
          - For config, it is the output of ``self.input_parser.explain_param(name)``.
          - For unknown, it is the ``Parameter`` from ``inspect.signature``.
        """

        #Mostly to avoid all possible spurious attribute errors
        getprovider = self.get_provider_func

        func = getprovider(provider_name)
        if isinstance(func, collect):
            def myfunc():
                pass
            name =  provider_name
            myfunc.__name__ = name
            innername = getattr(func.function, '__name__', str(func.function))
            myfunc.__doc__ = f"The result of `{innername}` for each in {func.fuzzyspec}."
            result = [myfunc]

            #TODO: This is incorrect. Fix it.
            if callable(func.function):
                func = func.function
            else:
                warnings.warn(f"Unimplemented help for {provider_name}")
                func = myfunc

        elif callable(func):
            result = [func]
        else:
            raise NotImplementedError(func)

        sig = inspect.signature(func)

        for param_name, param in sig.parameters.items():
            try:
                getprovider(param_name)
            except AttributeError:
                if self.input_parser.get_parse_func(param_name):
                    result.append(('config', param_name,
                        self.input_parser.explain_param(param_name)))
                elif self.input_parser.get_produce_func(param_name):
                    result.append(('produce', param_name,
                        self.input_parser.explain_param(param_name)))
                else:
                    result.append(('unknown', param_name, param))
            else:
                result.append(('provider', param_name,
                               self.explain_provider(param_name)))
        return result

    def expand_fuzzytarget_spec(self, fuzzytarget):
        name, fuzzy, nsroot ,extra_args = fuzzytarget
        specs = self.input_parser.process_fuzzyspec(fuzzy,
                                                self.rootns, parents=[name],
                                                initial_spec=nsroot)
        return specs


    def resolve_fuzzytargets(self):
        for target in self.fuzzytargets:
            self.resolve_fuzzytarget(target)

    def resolve_fuzzytarget(self, fuzzytarget):
        if not isinstance(fuzzytarget, FuzzyTarget):
            fuzzytarget = FuzzyTarget(*fuzzytarget)

        specs = self.expand_fuzzytarget_spec(fuzzytarget)

        for spec in specs:
            self.process_targetspec(fuzzytarget.name, spec, fuzzytarget.extraargs)

    def process_targetspec(self, name, nsspec, extraargs=None,
                            default=EMPTY):

        log.debug("Processing target %s" % name)

        gen = self._process_requirement(name, nsspec, extraargs=extraargs,
                                        default=default, parents=[])
        gen.send(None)
        try:
            gen.send(None)
        except StopIteration:
            pass
        else:
            raise RuntimeError()

    def _process_requirement(self, name, nsspec, *, extraargs=None,
                            default=EMPTY, parents=None):
        """Create nodes so as to satisfy the requirement specified by the
        arguments."""
        if parents is None:
            parents = []

        log.debug("Processing requirement: %s" % (name,))

        ns = namespaces.resolve(self.rootns, nsspec)
        if extraargs is None:
            extraargs = ()


        #First try to find the name in the namespace
        try:
            put_index, val = self.input_parser.resolve_key(name, ns, parents=parents, currspec=nsspec)
            log.debug("Found %s for spec %s at %s"%(name, nsspec, put_index))

        except InputNotFoundError as e:
            #See https://www.python.org/dev/peps/pep-3110/
            saved_exception = e
            #Handle this case later
            pass
        else:
            if extraargs:
                raise ResourceNotUnderstood(name, "The resource %s name is "
                "already present in the input, but some arguments were "
                "passed to compute it: %s" % (name, extraargs), parents[-1])

            if isinstance(val, ExplicitNode):
                yield from self._make_node((name, val.value), nsspec, extraargs, parents)
            else:
                yield put_index, val
            return

        #If the name is not in the providers, either it is an extra argument
        #or is missing

        if not self.is_provider_func(name):
            if default is EMPTY:
                raise saved_exception
            else:
                put_index = None
                yield put_index, default
                return

        #here we handle the case where the requirement is a provider and
        #make a new node for it.
        yield from self._make_node(name, nsspec, extraargs, parents)


    def _make_node(self, name_or_tuple, nsspec, extraargs, parents):
        """Make a node from the input arguments as well as any nodes required
        from that node."""
        if isinstance(name_or_tuple, tuple):
            name, f = name_or_tuple
        else:
            name = name_or_tuple
            f = self.get_provider_func(name)
        if isinstance(f, target_map):
            yield from self._make_collect_targets(f, name, nsspec, parents)
        elif isinstance(f, collect):
            yield from self._make_collect(f, name, nsspec, parents)
        else:
            yield from self._make_callspec(f, name, nsspec, extraargs, parents)

    def _create_default_key(self, name, nsspec, put_index=None, defaults=None):
        """Push a namespace level for a node to store input values.
        Return the nsspec of that node."""
        if defaults is None:
            defaults = {}
        if put_index is None:
            put_index = 0

        defaults_label = '_' + name + '_defaults'
        nsspec = (*nsspec[:len(nsspec)-put_index], defaults_label)
        parent_ns = namespaces.resolve(self.rootns, nsspec[:-1])
        namespaces.push_nslevel(parent_ns, defaults_label, defaults)
        return nsspec


    def _make_callspec(self, f, name, nsspec, extraargs, parents):
        """Make a normal node that calls a function."""

        defaults = {}
        s = inspect.signature(f)
        if(extraargs):
            defaults.update(dict(extraargs))

        #Note that this is the latest possible put_index and not len - 1
        #because there is also the root namespace.
        put_index = len(nsspec)
        gens = []
        for param_name, param in s.parameters.items():
            default = defaults.get(param_name, param.default)
            gen = self._process_requirement(param_name, nsspec, extraargs=None,
                                     default=default, parents=[name, *parents])
            index, _ = gen.send(None)
            log.debug("put_index for %s is %s" % (param_name, index))
            if index is None:
                defaults[param_name] = default
            elif index < put_index:
                put_index = index
            gens.append(gen)

        #The namespace stack (put_index) goes in the opposite direction
        #of the nsspec. put_index==len(nsspec)==len(ns.maps)-1
        #corresponds to the root namespace, and put_index=0 to the current
        #spec.

        #We need the len bit for the case put_index==0
        newnsspec = self._create_default_key(name, nsspec, put_index, defaults)
        log.debug("New spec for %s is: %s" %(name, newnsspec,))


        ns = namespaces.resolve(self.rootns, newnsspec)

        cs = CallSpec(f, tuple(s.parameters.keys()), name,
                      newnsspec)
        already_exists = cs in self.graph
        if already_exists:
            log.debug("Node '%s' already in the graph.", cs)
        else:
            log.debug("Appending node '%s'." % (cs,))
            self.graph.add_or_update_node(cs)
        for gen in gens:
            try:
               gen.send(cs)
            except StopIteration:
                pass
            else:
                raise RuntimeError()


        required_by = yield put_index, cs
        if required_by is None:
            outputs = set()
        else:
            outputs = set([required_by])
        self.graph.add_or_update_node(cs, outputs=outputs)

        #Do not repeat the checks for the same node
        if already_exists:
            return

        try:
            check_types(f, ns)
        except BadInputType as e:
            raise ResourceError(name, e, parents) from e

        if hasattr(f, 'checks'):
            for check in f.checks:
                try:
                    check(callspec=cs, ns=ns, graph=self.graph,
                          environment=self.environment)
                except CheckError as e:
                    raise ResourceError(name, e, parents) from e

    def _make_collect(self, f, name, nsspec, parents):
        """Make a node that spans a function over the values in a list and
        collects them in another list."""
        newparents = [name, *parents]

        myspec = self._create_default_key(name, nsspec)



        specs = self.input_parser.process_fuzzyspec(f.fuzzyspec,
                                                    self.rootns,
                                                    parents=newparents,
                                                    initial_spec=nsspec)
        myns = namespaces.resolve(self.rootns, myspec)
        myns[collect.resultkey] = OrderedDict.fromkeys(range(len(specs)))

        if not isinstance(f.function, str):
            newname = f.function.__name__
        else:
            newname = f.function

        compiletime = True

        collspec = CollectSpec(f, (), name, myspec)

        gens = []
        for i, spec in enumerate(specs):

            gen = self._process_requirement(newname, spec,
                                            parents=newparents)


            try:
                _, newcs = gen.send(None)
            except InputNotFoundError:
                if f.element_default is EMPTY:
                    raise
                newcs = f.element_default
            else:
                gens.append(gen)

            if isinstance(newcs, Node):
                flagargs = (('target', collspec), ('index', i))
                self._node_flags[newcs].add((add_to_dict_flag, flagargs))
                compiletime = False
            else:
                #TODO: Find a better way to do this: E.g. input nodes
                myns[collect.resultkey][i] = newcs

        if compiletime:
           value = f(self.rootns, myspec)
           namespaces.resolve(self.rootns, nsspec)[name] = value
           yield 0, value
        else:
            required_by = yield 0, collspec
            log.debug("Appending node {}".format(collspec))
            if required_by is None:
                outputs = set()
            else:
                outputs = set([required_by])
            value = collspec
            self.graph.add_or_update_node(collspec, outputs=outputs)
        for gen in gens:
            try:
                gen.send(value)
            except StopIteration:
                pass
            else:
                raise RuntimeError()

    def _make_collect_targets(self, colltargets, name, nsspec, parents):

        newparents = [name, *parents]

        myspec = self._create_default_key(name, nsspec)

        my_node = CollectMapSpec(colltargets, () ,name, myspec)

        required_by = yield 0, my_node

        if required_by is None:
            outputs = set()
        else:
            outputs = set([required_by])

        self.graph.add_or_update_node(my_node, outputs=outputs)


        def walk(d, spec):
            for target in d['targets']:
                target = target._replace(rootspec=spec)

                target_specs = self.expand_fuzzytarget_spec(target)
                for i, tspec in enumerate(target_specs):
                    gen = self._process_requirement(name=target.name, nsspec=tspec,
                                                    extraargs=target.extraargs, parents=newparents)
                    index, tnode = gen.send(None)
                    try:
                        gen.send(my_node)
                    except StopIteration:
                        pass
                    else:
                        raise RuntimeError()

            for fuzzy, newd in d['withs'].items():
                newspecs = self.input_parser.process_fuzzyspec(fuzzy,
                                                self.rootns, parents=newparents,
                                                initial_spec=spec)
                for newspec in newspecs:
                    walk(newd, newspec)

        walk(colltargets.root, nsspec)
