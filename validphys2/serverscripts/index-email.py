from collections import defaultdict
import json
import pickle
from urllib.parse import urljoin, urlparse

from bs4 import BeautifulSoup
import requests

ARCHIVES_URL = 'https://lists.cam.ac.uk/mailman/private/ucam-nnpdf/'

USERNAME = 'bot@zigzah.com'
PASSWORD_FILE = 'EMAIL_BOT_PASSWORD'


def make_login():
    with open(PASSWORD_FILE) as f:
        password = f.read().strip()
    return {"password": password, "username": USERNAME, "name": "submit"}


def make_soup(data):
    return BeautifulSoup(data, features="html5lib")


def get_archive_index(session):
    resp = session.post(ARCHIVES_URL, data=make_login())
    resp.raise_for_status()
    return resp.text


def get_thread_index(month_url, session):
    resp = session.get(month_url)
    resp.raise_for_status()
    return resp.text


def parse_threads(archive_index):
    soup = make_soup(archive_index)
    return [
        urljoin(ARCHIVES_URL, th.attrs['href']) for th in soup.find_all('a', string='[ Thread ]')
    ]


def parse_emails(thread_index, month_url):
    soup = make_soup(thread_index)
    return [
        urljoin(month_url, em.attrs['href'])
        for em in soup.find_all('a', attrs={'name': True, 'href': True})
    ]


def get_email(email_url, session):
    resp = session.get(email_url)
    resp.raise_for_status()
    return resp.text


def parse_email(email, email_url):
    res = {}

    def good_link(tag):
        if tag.name != 'a':
            return False
        if not tag.attrs['href'].startswith('https://vp.nnpdf.science/'):
            return False
        if any(p.name == 'blockquote' for p in tag.parents):
            return False
        return True

    soup = make_soup(email)
    links = soup.body.find_all(good_link, recursive=True)
    for link in links:
        p = urlparse(link.attrs['href']).path
        fragments = p.split('/')
        if len(fragments) >= 2:
            res[fragments[1]] = [email_url, str(soup.title.string)]
    return res


if __name__ == '__main__':
    try:
        with open('email_mentions.json') as f:
            res = defaultdict(list, json.load(f))
    except FileNotFoundError:
        res = defaultdict(list)

    try:
        with open('seen_emails_cache.pkl', 'rb') as f:
            seen_data = pickle.load(f)
        seen_months = seen_data['seen_months']
        seen_emails = seen_data['seen_emails']
    except FileNotFoundError:
        seen_months = set()
        seen_emails = set()

    s = requests.Session()
    idx = get_archive_index(session=s)
    for mindex, month_url in enumerate(parse_threads(idx)):
        if month_url in seen_months:
            continue
        # Could still add emails to last month
        if mindex != 0:
            seen_months.add(month_url)
        thindex = get_thread_index(month_url, s)
        for email_url in parse_emails(thindex, month_url):
            if email_url in seen_emails:
                continue
            seen_emails.add(email_url)
            email = get_email(email_url, s)
            email_res = parse_email(email, email_url)
            for k, v in email_res.items():
                res[k].append(v)
    with open('email_mentions.json', 'w') as f:
        json.dump(res, f)
    with open('seen_emails_cache.pkl', 'wb') as f:
        pickle.dump({'seen_months': seen_months, 'seen_emails': seen_emails}, f)
