# -*- coding: utf-8 -*-
"""
Generate an index with the existing internal or unpublished PDFS.
"""

import json
import pathlib

root = '/home/nnpdf/WEB/pdfs'

glob = '*gz'

indexname = 'pdfdata.json'

if __name__ == '__main__':
    p = pathlib.Path(root)
    files = p.glob(glob)
    files = [f.name for f in files]
    with (p / indexname).open('w') as f:
        json.dump({'files': files}, f, separators=(',', ':'))
