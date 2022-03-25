# -*- coding: utf-8 -*-
"""
Generate an index with the existing fits
"""

import pathlib
import json
import itertools
import tarfile
import sys
import pickle

import ruamel_yaml as yaml

root = 'WEB/fits'
seen_cache_path = "seen_fits.pkl"

#TODO: Find a better way
glob1 = '*.tar.gz'
glob2 = '*.tgz'

indexname = 'fitdata.json'

def true_stem(p):
    # Note: See Path.prefixes to cross check implementation.
    return p.name.lstrip(".").split(".", 1)[0]

def process_tar_meta(p):
    with tarfile.open(p) as tf:
        ef = tf.extractfile(f"{true_stem(p)}/filter.yml")
        meta = yaml.YAML().load(ef)
    return meta

def process_file_iter(files, seen_cache):
    names = []
    metas = {}
    for p in files:
        names.append(p.name)
        fitname = true_stem(p)
        if p.name in seen_cache:
            meta = seen_cache[fitname]
        else:
            try:
                meta = process_tar_meta(p)
            except Exception as e:
                print(f"Failed processing {p}", e, file=sys.stderr)
            else:
                metas[fitname] = meta
    return names, metas


if __name__ == '__main__':
    p = pathlib.Path(root)
    files = itertools.chain(p.glob(glob1), p.glob(glob2))

    try:
        with open(seen_cache_path, "rb") as f:
            seen_cache = pickle.load(f)
    except FileNotFoundError as e:
        seen_cache = {}


    names, metas = process_file_iter(files, seen_cache)

    with (p/indexname).open('w') as f:
        json.dump({'files':files, 'meta': metas}, f, indent=4)

    with open(seen_cache_path, "wb") as f:
        pickle.dump(metas, f)
