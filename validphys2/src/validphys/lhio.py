# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 22:00:01 2015

@author: zah

A module that reads and writes LHAPDF grids.
"""

import os
import os.path as osp
import sys
import shutil

import numpy as np
import pandas as pd
import yaml

import applwrap

from smpdflib import lhaindex

def split_sep(f):
    for line in f:
        if line.startswith(b'---'):
            break
        yield line

def read_xqf_from_file(f):

    lines = split_sep(f)
    try:
        (xtext, qtext, ftext) = [next(lines) for _ in range(3)]
    except StopIteration:
        return None
    xvals = np.fromstring(xtext, sep = " ")
    qvals = np.fromstring(qtext, sep = " ")
    fvals = np.fromstring(ftext, sep = " ", dtype=np.int)
    vals = np.fromstring(b''.join(lines), sep= " ")
    return pd.Series(vals, index = pd.MultiIndex.from_product((xvals, qvals, fvals)))


def read_xqf_from_lhapdf(pdf, replica, rep0grids):
    indexes = tuple(rep0grids.index)
    vals = []
    for x in indexes:
        vals += [pdf.xfxQ(replica,x[3],x[1],x[2])]
    return pd.Series(vals, index = rep0grids.index)

def read_all_xqf(f):
    while True:
        result = read_xqf_from_file(f)
        if result is None:
            return
        yield result

def load_replica( pdf, rep, rep0grids=None):

    sys.stdout.write("-> Reading replica from LHAPDF %d \r" % rep)
    sys.stdout.flush()

    suffix = str(rep).zfill(4)

    pdf_name = str(pdf)

    path = osp.join(lhaindex.finddir(pdf_name),
                    pdf_name + "_" + suffix + ".dat")

    with open(path, 'rb') as inn:
        header = b"".join(split_sep(inn))

        if rep0grids is not None:
            xfqs = read_xqf_from_lhapdf(pdf, rep, rep0grids)
        else:
            xfqs = list(read_all_xqf(inn))
            xfqs = pd.concat(xfqs, keys=range(len(xfqs)))
    return header, xfqs

#Split this to debug easily
def _rep_to_buffer(out, header, subgrids):
    sep = b'---'
    out.write(header)
    out.write(sep)
    for _,g in subgrids.groupby(level=0):
        out.write(b'\n')
        ind = g.index.get_level_values(1).unique()
        np.savetxt(out, ind, fmt='%.7E',delimiter=' ', newline=' ')
        out.write(b'\n')
        ind = g.index.get_level_values(2).unique()
        np.savetxt(out, ind, fmt='%.7E',delimiter=' ', newline=' ')
        out.write(b'\n')
        #Integer format
        ind = g.index.get_level_values(3).unique()
        np.savetxt(out, ind, delimiter=' ', fmt="%d",
                      newline=' ')
        out.write(b'\n ')
        #Reshape so printing is easy
        reshaped = g.reshape((len(g.groupby(level=1))*len(g.groupby(level=2)),
                              len(g.groupby(level=3))))
        np.savetxt(out, reshaped, delimiter=" ", newline="\n ", fmt='%14.7E')
        out.write(sep)

def write_replica(rep, pdf_name, header, subgrids):
    suffix = str(rep).zfill(4)
    with open(pdf_name + "_" + suffix + ".dat", 'wb') as out:
        _rep_to_buffer(out, header, subgrids)

def load_all_replicas(pdf, db=None):
    if db is not None:
        #removing str() will crash as it casts to unicode due to pdf name
        key = str("(load_all_replicas, %s)" % pdf.get_key())
        if key in db:
            return db[key]
    rep0headers, rep0grids = load_replica(pdf, 0)

    headers, grids = zip(*[load_replica(pdf, rep, rep0grids)
                         for rep in range(1, len(pdf))])
    result = [rep0headers] + list(headers), [rep0grids] + list(grids)
    if db is not None:
        db[key] = result
    return result

def big_matrix(gridlist):
    central_value = gridlist[0]
    X = pd.concat(gridlist[1:], axis=1,
                 keys=range(1,len(gridlist)+1), #avoid confusion with rep0
                 ).subtract(central_value, axis=0)
    if np.any(X.isnull()) or X.shape[0] != len(central_value):
        raise ValueError("Incompatible grid specifications")
    return X

def hessian_from_lincomb(pdf, V, set_name=None, folder = None, db=None,
                         extra_fields=None):
    """Construct a new LHAPDF grid from a linear combination of members"""

    # preparing output folder
    neig = V.shape[1]

    base = applwrap.getlhapdfpath()[-1] + "/" + str(pdf) + "/" + str(pdf)
    if set_name is None:
        set_name = str(pdf) + "_hessian_" + str(neig)
    if folder is None:
        folder = ''
    set_root = os.path.join(folder,set_name)
    if not os.path.exists(set_root): os.makedirs(os.path.join(set_root))

    # copy replica 0
    shutil.copy(base + "_0000.dat", set_root + "/" + set_name + "_0000.dat")

    with open(base + ".info", 'r') as inn, \
         open(set_root + "/" + set_name + ".info", 'w') as out:

        for l in inn.readlines():
            if l.find("SetDesc:") >= 0:
                out.write("SetDesc: \"Hessian " + str(pdf) + "_hessian\"\n")
            elif l.find("NumMembers:") >= 0:
                out.write("NumMembers: " + str(neig+1) + "\n")
            elif l.find("ErrorType: replicas") >= 0:
                out.write("ErrorType: symmhessian\n")
            else:
                out.write(l)
        if extra_fields is not None:
            yaml.dump(extra_fields, out, default_flow_style=False)

    headers, grids = load_all_replicas(pdf, db=db)
    hess_name = set_root + '/' + set_name
    result  = (big_matrix(grids).dot(V)).add(grids[0], axis=0, )
    hess_header = b"PdfType: error\nFormat: lhagrid1\n"
    for column in result.columns:
        write_replica(column + 1, hess_name, hess_header, result[column])

    return set_root
