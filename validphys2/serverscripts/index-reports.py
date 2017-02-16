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

    date = datetime.date.fromtimestamp(p.stat().st_mtime).isoformat()
    if not title:
        title = "Validphys output (untitled)"

    titlelink = '<a href="%s">%s</a>' % (url, title)
    return (titlelink, author, date, tags)

def make_index(root_path, out):
    root_path = pathlib.Path(root_path)
    data = [register(p) for p in root_path.iterdir() if p.is_dir()]
    with open(out, 'w') as f:
        json.dump({'data':data, }, f)

def main():
    make_index(ROOT, OUT)

if __name__ == '__main__':
    main()


