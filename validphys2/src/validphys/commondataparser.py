"""
This module implements parsers for vommondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
import io
import functools
import tarfile
import dataclasses

import numpy as np
import pandas as pd

from validphys.coredata import CommonData 
