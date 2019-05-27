"""
API for resource engine
"""

import logging
import importlib

from reportengine.resourcebuilder import ResourceBuilder
from reportengine.resourcebuilder import FuzzyTarget

log = logging.getLogger(__name__)

class API:
    """The API class"""
    def __init__(self, providers, config_cls, env_cls, **kwargs):
        prov_list = []
        for prov in providers:
            if isinstance(prov, str):
                try:
                    mod = importlib.import_module(prov)
                except ImportError:
                    log.error("Could not import module %s", prov)
                    raise
            else:
                mod = prov
            prov_list.append(mod)

        self.provider_loaded = prov_list
        self.config_class = config_cls
        self.loadedenv = env_cls(**kwargs)

    def __call__(self, actions: str, **kwargs):
        fuzzytarg = [FuzzyTarget(actions, (), (), ())]
        c = self.config_class(kwargs, environment=self.loadedenv)
        builder = ResourceBuilder(c, self.provider_loaded, fuzzytarg, perform_final=False)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        res = builder.rootns[actions]
        return res

    def __getattr__(self, name):
        def closure(**kwargs):
            return self.__call__(name, **kwargs)
        return closure
