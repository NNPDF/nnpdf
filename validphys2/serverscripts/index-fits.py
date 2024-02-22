# -*- coding: utf-8 -*-
"""
Generate an index with the existing fits
"""

import itertools
import json
import pathlib

root = '/home/nnpdf/WEB/fits'

# TODO: Find a better way
glob1 = '*.tar.gz'
glob2 = '*.tgz'

indexname = 'fitdata.json'

if __name__ == '__main__':
    p = pathlib.Path(root)
    files = itertools.chain(p.glob(glob1), p.glob(glob2))

    files = [f.name for f in files]
    with (p / indexname).open('w') as f:
        json.dump({'files': files}, f, indent=4)
