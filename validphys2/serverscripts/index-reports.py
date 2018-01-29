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
from collections import ChainMap, defaultdict

import ruamel_yaml as yaml
from bs4 import BeautifulSoup

ROOT = '/home/nnpdf/WEB/validphys-reports'
ROOT_URL = 'http://pcteserver.mi.infn.it/~nnpdf/validphys-reports/'
OUT = '/home/nnpdf/WEB/validphys-reports/index.json'

EMPTY = '-'

DEFAULTS = dict(author=EMPTY, title=EMPTY, keywords=[])


REQUIRED_FILE_METADATA = {'title', 'author', 'keywords'}

def meta_from_html(f):
    soup = BeautifulSoup(f, 'lxml')
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
    #As it turns out, 'soup.title.string' idiotically doesn't
    #return a strig but rather an object with the reference to
    #the whole parse tree, causing a huge memory leak.
    return dict(title=str(title), author=author, keywords=tags)

class TagProps():
    def __init__(self, count=0, last_timestamp=0):
        self.count = count
        self.last_timestamp = last_timestamp

    __slots__ = ('count', 'last_timestamp')

def meta_from_path(p):
    meta = ChainMap(DEFAULTS)
    yaml_meta = p/'meta.yaml'
    if yaml_meta.exists():
        with yaml_meta.open() as f:
            yaml_res = yaml.safe_load(f)
    else:
        yaml_res = {}
    index = p/'index.html'
    #Only do the expensive HTML parsing if we actually need a key
    if REQUIRED_FILE_METADATA - yaml_res.keys() and index.exists():
        with index.open() as f:
            meta = meta.new_child(meta_from_html(f))
    meta = meta.new_child(yaml_res)
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
    keywords = defaultdict(TagProps)
    for p in root_path.iterdir():
        if p.is_dir():
            res = register(p)
            data.append(res)
            newkeywords = res[3]
            timestamp = res[2][1]
            for k in newkeywords:
                props = keywords[k]
                props.count+=1
                props.last_timestamp = max(props.last_timestamp, timestamp)

    keylist = sorted(keywords.items(), key=lambda x: -x[1].last_timestamp)
    keywordmap = [(k, v.count) for k,v in keylist]


    with open(out, 'w') as f:
        json.dump({'data':data, 'keywords':keywordmap}, f)

def main():
    make_index(ROOT, OUT)

if __name__ == '__main__':
    main()


