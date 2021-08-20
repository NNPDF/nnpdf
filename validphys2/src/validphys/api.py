"""
api.py

This module contains the `reportengine` programmatic API, initialized with the
validphys providers, Config and Environment.

Example:
--------

Simple Usage:

>> from validphys.api import API
>> fig = API.plot_pdfs(pdf="NNPDF_nlo_as_0118", Q=100)
>> fig.show()

"""
import logging

from reportengine import api

from validphys.app import providers
from validphys.config import Config, Environment

log = logging.getLogger(__name__)

# API needed its own module, so that it can be used with any Matplotlib backend
# without breaking validphys.app
API = api.API(providers, Config, Environment)
