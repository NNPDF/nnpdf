"""
api.py

Interface module which allows the user to access the necessary reportengine technology
to call validphys providers with the same inputs as a validphys runcard from a python interpreter.

Additionally the *provider modules* that serve as source to the validphys actions are
declared here.

Standard use case is something like:

from validphys.api import API

vp2 = API()
figure = vp2.plot_pdfs({'pdfs':['NNPDF31_nlo_as_0118], 'Q':1.65})


"""
import logging
import importlib

from validphys.config import Environment, Config
from reportengine.resourcebuilder import ResourceBuilder
from reportengine.resourcebuilder import FuzzyTarget

log = logging.getLogger(__name__)

providers = [
    'validphys.results',
    'validphys.pdfgrids',
    'validphys.pdfplots',
    'validphys.dataplots',
    'validphys.fitdata',
    'validphys.arclength',
    'validphys.sumrules',
    'validphys.reweighting',
    'validphys.kinematics',
    'validphys.correlations',
    'validphys.chi2grids',
    'validphys.eff_exponents',
    'validphys.paramfits.dataops',
    'validphys.paramfits.plots',
    'validphys.theorycovariance',
    'validphys.replica_selector',
    'validphys.MCgen_checks',
    'validphys.closure',
    'validphys.api',
    'reportengine.report']

class API:
    config_class = Config
    environment_class = Environment

    def __init__(self, **kwargs):
        prov_list = []
        for prov in providers:
            try:
                mod = importlib.import_module(prov)
            except ImportError:
                log.error("Could not import module %s", mod)
                raise
            prov_list.append(mod)
        self.provider_loaded = prov_list
        self.loaderenv = Environment(**kwargs)

    def __call__(self, actions: str, **kwargs):
        fuzzytarg = [FuzzyTarget(actions, (), (), ())]
        c = self.config_class(kwargs, environment=self.loaderenv)
        builder = ResourceBuilder(c, self.provider_loaded, fuzzytarg)
        builder.resolve_fuzzytargets()
        builder.execute_sequential(perform_final=False)
        res = builder.rootns[actions]
        return res

    def __getattr__(self, name):
        print(type(name))
        def closure(**kwargs):
            return self.__call__(name, **kwargs)
        return closure
