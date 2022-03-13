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
import sys
import traceback
from collections import ChainMap, defaultdict

import ruamel_yaml as yaml
from bs4 import BeautifulSoup
#TODO: Move the thumbnail logic somewhere
import skimage.transform
import skimage.io
import numpy as np

ROOT = '/home/nnpdf/validphys-reports'
ROOT_URL = 'https://vp.nnpdf.science/'
OUT = '/home/nnpdf/validphys-reports/index.json'
THUMBNAILS = '/home/nnpdf/thumbnails/'
EMAIL_MENTIONS_FILE = 'email_mentions.json'

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
    #'soup.title.string' doesn't
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
    yaml_res = {}
    if yaml_meta.exists():
        with yaml_meta.open() as f:
            try:
                yaml_res = yaml.safe_load(f)
            except yaml.YAMLError as e:
                print(f"Error processing {yaml_meta}: {e}", file=sys.stderr)
    index = p/'index.html'
    #Only do the expensive HTML parsing if we actually need a key
    if REQUIRED_FILE_METADATA - yaml_res.keys() and index.exists():
        with index.open() as f:
            meta = meta.new_child(meta_from_html(f))
    meta = meta.new_child(yaml_res)
    return meta

def make_single_thumbnail(f, shape=(100, 150)):
    img = skimage.io.imread(f)
    res = skimage.transform.resize(
        img, shape, anti_aliasing=True, mode='constant')
    return res

def make_4_img_thumbnail(paths, shape=(100, 150)):
    w, h = shape
    whalf, hhalf = w // 2, h // 2
    positions = (
        (slice(0,whalf), slice(0,hhalf)),
        (slice(whalf,w), slice(0,hhalf)),
        (slice(0,whalf), slice(hhalf,h)),
        (slice(whalf,w), slice(hhalf,h))
    )
    res = np.zeros((*shape, 4))
    imgs = skimage.io.imread_collection(paths)
    for img, pos in zip(imgs, positions):
        res[pos] = skimage.transform.resize(
            img, (whalf, hhalf), anti_aliasing=True, mode='constant')
    return res

def make_thumbnail(folder):
    folder = pathlib.Path(folder)
    pngs = sorted(folder.glob('*.png'))
    if not pngs:
        return None
    if len(pngs) < 4:
        return make_single_thumbnail(pngs[0])
    else:
        l = len(pngs)
        imgs = pngs[:l-(l%4):l//4]
        return make_4_img_thumbnail(imgs)


def thumbnail_tag(name):
    return f'{ROOT_URL}thumbnails/{name}"'

def handle_thumbnail(p):
    dest = (pathlib.Path(THUMBNAILS) / p.name).with_suffix('.png')
    name = dest.name
    if dest.exists():
        return thumbnail_tag(name)
    figures = (p / 'figures')
    if figures.is_dir():
        try:
            res = make_thumbnail(figures)
            if res is not None:
                skimage.io.imsave(dest, res)
                return thumbnail_tag(name)
        except Exception as e:
            print("Could not process thumbnails in", figures, e, file=sys.stderr)
            traceback.print_exc()
            return None
    return None

def register(p, emails):
    path_meta = meta_from_path(p)
    title, author, tags = path_meta['title'], path_meta['author'], path_meta['keywords']
    url = ROOT_URL + p.name

    #Use the timestamp for sorting and the string for displaying
    timestamp = p.stat().st_mtime
    date = datetime.date.fromtimestamp(timestamp).isoformat()
    if not title or not isinstance(title, str):
        title = "Validphys output (untitled)"

    if isinstance(tags, str):
        tags = [tags]
    if not isinstance(tags, list):
        tags = []
    if not isinstance(author, str):
        author = "<INVALID AUTHOR>"

    emaillinks = ' '.join(
        f'<a href="{url}", title="{title}">📧</a>' for (url, title) in emails
    )

    titlelink = f'<a href="{url}">{title}</a> {emaillinks}'

    thumbnail = handle_thumbnail(p)

    return (titlelink, author, [date, timestamp], tags, thumbnail)


def get_all_emails():
    try:
        with open(EMAIL_MENTIONS_FILE) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def make_index():
    root_path = pathlib.Path(ROOT)
    emails = get_all_emails()
    data = []
    keywords = defaultdict(TagProps)
    for p in root_path.iterdir():
        if p.is_dir():
            try:
                res = register(p, emails.get(p.name, []))
                data.append(res)
                newkeywords = res[3]
                timestamp = res[2][1]
                for k in newkeywords:
                    props = keywords[k]
                    props.count+=1
                    props.last_timestamp = max(props.last_timestamp, timestamp)
            except:
                print("Error processing folder", p,file=sys.stderr)
                raise

    keylist = sorted(keywords.items(), key=lambda x: -x[1].last_timestamp)
    keywordmap = [(k, v.count) for k,v in keylist]


    with open(OUT, 'w') as f:
        json.dump({'data':data, 'keywords':keywordmap}, f)


if __name__ == '__main__':
    make_index()
