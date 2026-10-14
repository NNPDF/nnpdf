"""INSPIRE API client to fetch bibtex entries."""

import dataclasses
from functools import lru_cache
import logging
import re
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

INSPIRE_API = "https://inspirehep.net/api/literature"
TIMEOUT = 10
MAX_RETRIES = 2
# Delay set to 5s to ensure it resets the inspire API limit
RETRY_DELAY = 5.0

HEADERS = {"User-Agent": "nnpdf-data-cli/1.0"}

_EXTRACTOR = re.compile(r"@\w+\{([^,]+),")


class InspireError(Exception):
    pass


@dataclasses.dataclass(frozen=True)
class BibliographyEntry:
    """A parsed BibTeX bibliography entry."""

    key: str
    raw_bibtex: str


def _extract_key(raw: str) -> str | None:
    """Extract BibTeX key from raw entry."""
    m = _EXTRACTOR.match(raw.strip())
    return m.group(1).strip() if m else None


@lru_cache
def _fetch_bibtex(record_id: int):
    """Fetch bibtex using the inspire API given a record id."""
    url = f"{INSPIRE_API}/{record_id}?format=bibtex"
    for attempt in range(MAX_RETRIES + 1):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=TIMEOUT) as resp:
                raw = resp.read().decode("utf-8")
            key = _extract_key(raw)
            if key is not None:
                logger.info(f"Received data for {key}: {record_id}")
                return BibliographyEntry(key, raw)
        except HTTPError as e:
            if e.code != 429 and e.code < 500:
                # We got a bona-fide error, probably something's wrong
                # unless it is 429, then we simply got rate-limited, just waiting should reset it
                logger.warning("iNSPIRE API returned HTTP %s for record %s", e.code, record_id)
                raise InspireError
        except URLError as e:
            # Most likely our fault, wait for a bit then retry
            pass
        time.sleep(RETRY_DELAY)
    raise InspireError


def fetch_bibtex(record_int: int):
    """Wrapper around _fetch_bibtex to convert known exceptions into Nones"""
    try:
        return _fetch_bibtex(record_int)
    except InspireError:
        return None
