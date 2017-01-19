import pathlib
import datetime
import json

from bs4 import BeautifulSoup

ROOT = '/home/nnpdf/WEB/validphys-reports'
ROOT_URL = 'http://pcteserver.mi.infn.it/~nnpdf/validphys-reports/'
OUT = '/home/nnpdf/WEB/validphys-reports/index.json'

EMPTY = '-'

def meta_from_html(f):
    soup = BeautifulSoup(f, 'html.parser')
    try:
        title = soup.title.string
    except Exception:
        title = EMPTY
    if not title:
        title = EMPTY
    try:
        author = soup.find('meta', {'name':'author'})['content']
    except Exception:
        author = EMPTY
    return title, author

def register(p):
    index = p/'index.html'
    if index.exists():
        title, author = meta_from_html(index.open())
    else:
        title = author = EMPTY
    url = ROOT_URL + p.name

    urllink = '<a href="%s">%s</a>' % (url, p.name)
    titlelink = '<a href="%s">%s</a>' % (url, title)

    date = datetime.date.fromtimestamp(p.stat().st_mtime).isoformat()

    return (titlelink, author, date, urllink)

def make_index(root_path, out):
    root_path = pathlib.Path(root_path)
    data = [register(p) for p in root_path.iterdir() if p.is_dir()]
    with open(out, 'w') as f:
        json.dump({'data':data, }, f)

def main():
    make_index(ROOT, OUT)

if __name__ == '__main__':
    main()


