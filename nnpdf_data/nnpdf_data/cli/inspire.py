"""INSPIRE API client to fetch bibtex entries."""

import dataclasses
import logging
import re
import time
from urllib.error import URLError
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

INSPIRE_API = "https://inspirehep.net/api/literature"
TIMEOUT = 10
MAX_RETRIES = 2
RETRY_DELAY = 1.0

HEADERS = {"User-Agent": "nnpdf-data-cli/1.0"}

_EXTRACTOR = re.compile(r"@\w+\{([^,]+),")


@dataclasses.dataclass(frozen=True)
class BibliographyEntry:
    """A parsed BibTeX bibliography entry."""

    key: str
    raw_bibtex: str


def _extract_key(raw: str) -> str | None:
    """Extract BibTeX key from raw entry."""
    m = _EXTRACTOR.match(raw.strip())
    return m.group(1).strip() if m else None


def fetch_bibtex(record_id: int):
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
        except URLError as e:
            code = getattr(e, "code", None)
            if code is not None and code < 500:
                # We got a bona-fide error, probably something's wrong
                return None
        time.sleep(RETRY_DELAY * (attempt + 1))
    return None
