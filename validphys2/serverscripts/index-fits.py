# -*- coding: utf-8 -*-
"""
Generate an index with the existing fits
"""

import pathlib
import json

root = '/home/nnpdf/WEB/fits'

glob = '*.tgz'

indexname = 'fitdata.json'

if __name__ == '__main__':
    p = pathlib.Path(root)
    files = p.glob(glob)
    files = [f.name for f in files]
    with (p/indexname).open('w') as f:
        json.dump({'files':files}, f, indent=4)

