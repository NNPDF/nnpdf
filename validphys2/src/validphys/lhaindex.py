# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Fri Jan 23 12:11:23 2015

@author: zah
"""
import fnmatch
from functools import lru_cache
import glob
import os
import os.path as osp
import re

import lhapdf

from reportengine.compat import yaml

_indexes_to_names = None
_names_to_indexes = None


def expand_index_names(globstr):
    return fnmatch.filter(get_names_to_indexes().keys(), globstr)


def expand_local_names(globstr):
    paths = get_lha_paths()
    return [
        name
        for path in paths
        for name in glob.glob1(path, globstr)
        if osp.isdir(osp.join(path, name))
    ]


def expand_names(globstr):
    """Return names of installed PDFs. If none is found,
    return names from index"""
    names = expand_local_names(globstr)
    if not names:
        names = expand_index_names(globstr)
    return names


def get_indexes_to_names():
    global _indexes_to_names
    if _indexes_to_names is None:
        _indexes_to_names = parse_index(get_index_path())
    return _indexes_to_names


def finddir(name):
    for path in get_lha_paths():
        d = osp.join(path, name)
        if osp.isdir(d):
            return d
    raise FileNotFoundError(name)


def isinstalled(name):
    """Check that name exists in LHAPDF dir"""
    return name and any(osp.isdir(osp.join(path, name)) for path in get_lha_paths())


def get_names_to_indexes():
    global _names_to_indexes
    if _names_to_indexes is None:
        _names_to_indexes = {name: index for index, name in get_indexes_to_names().items()}
    return _names_to_indexes


def get_pdf_indexes(name):
    """Get index in the amc@nlo format"""
    info = parse_info(name)
    ind = info['SetIndex']
    num_members = info['NumMembers']
    return {
        'lhapdf_id': ind,
        'lhapdf_min': ind + (num_members > 1),
        'lhapdf_max': ind + num_members - 1,
    }


def get_pdf_name(index):
    return get_indexes_to_names()[str(index)]


def parse_index(index_file):
    d = {}
    name_re = '(\d+)\s+(\S+)'
    with open(index_file) as localfile:
        for line in localfile.readlines():
            m = re.match(name_re, line)
            index = m.group(1)
            d[index] = m.group(2)
    return d


def get_collaboration(name):
    try:
        col = name[: name.index('_')]
    except ValueError:
        col = name
    return col


def as_from_name(name):
    """Annoying function needed because this is not in the info files. as(M_z)
    there is actually as(M_ref)."""
    match = re.search(r'as[_]?([^_\W]+)\_?', name)
    if match:
        return match.group(1)
    else:
        return name


def infofilename(name):
    for path in get_lha_paths():
        info = osp.join(path, name, name + '.info')
        if osp.exists(info):
            return info
    raise FileNotFoundError(name + ".info")


@lru_cache()
def parse_info(name):
    with open(infofilename(name)) as infofile:
        result = yaml.safe_load(infofile)
    return result


def get_lha_paths():
    return lhapdf.paths()


def get_lha_datapath():
    return get_lha_paths()[-1]


def get_index_path(folder=None):
    if folder:
        index_file = os.path.join(folder, 'pdfsets.index')
    if folder is None or not osp.exists(index_file):
        folder = get_lha_datapath()
    index_file = os.path.join(folder, 'pdfsets.index')
    return index_file
