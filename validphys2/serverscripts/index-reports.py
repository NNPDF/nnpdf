"""index-reports.py
Read the designed report folder and generate a JSON index to be
consumed by the index table (and possibly other applications).

The fields are infered from the import file attributes (time), a file
called meta.yaml in the report folder and finally the html
attributes, in that order.
"""
import pathlib
import datetime
import json
import re
from collections import ChainMap

import ruamel_yaml as yaml
from bs4 import BeautifulSoup

ROOT = '/home/nnpdf/WEB/validphys-reports'
ROOT_URL = 'http://pcteserver.mi.infn.it/~nnpdf/validphys-reports/'
OUT = '/home/nnpdf/WEB/validphys-reports/index.json'

EMPTY = '-'

DEFAULTS = dict(author=EMPTY, title=EMPTY, keywords=[])

def meta_from_html(f):
    soup = BeautifulSoup(f, 'html.parser')
    try:
        title = soup.title.string
    except Exception:
        title = None
    try:
        author = soup.find('meta', {'name':'author'})['content']
    except Exception:
        author = EMPTY

    try:
        tagtext = soup.find('meta', {'name':'keywords'})['content']
    except Exception:
        tags = []
    else:
        tags = re.split(r"\s*,\s*", tagtext)
    return dict(title=title, author=author, keywords=tags)

def meta_from_path(p):
    meta = ChainMap(DEFAULTS)
    index = p/'index.html'
    if index.exists():
        with index.open() as f:
           meta = meta.new_child(meta_from_html(f))
    yaml_meta = p/'meta.yaml'
    if yaml_meta.exists():
        with yaml_meta.open() as f:
            meta = meta.new_child(yaml.safe_load(f))
    return meta


def register(p):
    path_meta = meta_from_path(p)
    title, author, tags = path_meta['title'], path_meta['author'], path_meta['keywords']
    url = ROOT_URL + p.name

    #Use the timestamp for sorting and the string for displaying
    timestamp = p.stat().st_mtime
    date = datetime.date.fromtimestamp(timestamp).isoformat()
    if not title:
        title = "Validphys output (untitled)"

    titlelink = '<a href="%s">%s</a>' % (url, title)
    return (titlelink, author, [date, timestamp], tags)

def make_index(root_path, out):
    root_path = pathlib.Path(root_path)
    data = []
    keywords = set()
    for p in root_path.iterdir():
        if p.is_dir():
            res = register(p)
            data.append(res)
            keywords.update(res[3])

    data = [register(p) for p in root_path.iterdir() if p.is_dir()]
    with open(out, 'w') as f:
        json.dump({'data':data, 'keywords':list(keywords)}, f)

def main():
    make_index(ROOT, OUT)

if __name__ == '__main__':
    main()


