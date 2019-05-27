# -*- coding: utf-8 -*-
"""
Nice formatting strings from namespace elements, to be used as path names.

Created on Mon May 16 12:48:20 2016

@author: Zahari Kassabov
"""
import re
import logging
from collections.abc import Collection

from reportengine import namespaces

log = logging.getLogger(__name__)

def normalize_name(name):
    """Remove characters that don't go well in a filename.
    That is, everything non alphanumerical, except from '_' and '-'."""
    return re.sub(r'[^\w_-]', '', str(name))


def get_nice_name(ns, nsspec, suffix=None):
    """Get a name by quering the parts of a namespace specification.
    ``ns`` should be a namespace ChainMap and ``nsspec`` a
    tuple with a valid specification
    (see the ``namespaces`` documentation for more details)"""
    parts = []
    currspec = []
    currns = ns
    for ele in nsspec:
        currspec.append(ele)
        val = namespaces.value_from_spcec_ele(currns, ele)


        #kind of ugly, but we don't want to dumpt compound types, and too long
        #filenames)
        if isinstance(val, Collection):
            val = str(ele)
        else:
            try:
                val = str(val)
            except Exception as e:
                log.debug("Could not convert a value (%r) to string: %s" %
                      (val, e))
                val = str(ele)
            else:
                if len(val) > 25:
                    val = str(ele)

        parts.append(normalize_name(val))

        currns = namespaces.resolve(ns, currspec)

    if suffix:
        parts.append(suffix)

    return '_'.join(parts)

def spec_to_nice_name(ns, nodespec, suffix=''):
    """Get a nice name from a ResourceExecutor node spec"""

    name_suffix = str(nodespec.function.__name__)
    if suffix:
        name_suffix +=  '_' + suffix

    #Remove the internal default namespace
    name = get_nice_name(ns, nodespec.nsspec[:-1], name_suffix)
    return name

