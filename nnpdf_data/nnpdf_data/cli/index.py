"""
Index the whole set of data plus a few nice utilities.
"""

import dataclasses
import fnmatch
import logging
from pathlib import Path

from validobj import ValidationError

from ..commondataparser import parse_set_metadata
from ..utils import DEFAULT_PATH_VPDATA, get_nnpdf_profile

logger = logging.getLogger(__name__)


def _discover_metadata_files() -> list[Path]:
    """Loop over all data files available in all paths declared in the
    NNPDF ``data_path`` profile. If a dataset appears more than once,
    only the one that appears first (i.e., with greater priority) is read.
    """
    seen = []
    ret = []
    for data_root in get_nnpdf_profile().get("data_path"):
        for metadata_file in data_root.glob("*/metadata.yaml"):
            if (setname := metadata_file.parent.name) in seen:
                continue
            seen.append(setname)
            ret.append(metadata_file)
    return ret


def build_index(skip_failures=True, skip_special=True):
    """Build the dataset index by scanning and loading all possible metadata files
    found in all possible data paths.

    If any metadata file cannot be read, simply skip it.
    """
    files = _discover_metadata_files()
    all_datasets = {}
    for metadata_file in files:
        try:
            set_metadata = parse_set_metadata(metadata_file)
            for obs in set_metadata.allowed_observables:
                observable = set_metadata.select_observable(obs)
                if observable.is_nnpdf_special and skip_special:
                    break
                if (dname := observable.name) in all_datasets:
                    raise ValueError(f"A collission should not be possible! {dname}")
                all_datasets[dname] = observable
        except ValidationError as e:
            logger.warning(f"Failed to parse {metadata_file} (error: {e})")
            if not skip_failures:
                raise e

    return dict(sorted(all_datasets.items()))


def _expand_pattern(pattern: str):
    """Expand brace expressions à la bash -> {a,b}* == [a*, b*]."""
    prefix, sep, rest = pattern.partition("{")
    if not sep:
        return [pattern]
    pattern_options, sep, suffix = rest.partition("}")
    if not sep:
        return [pattern]
    return [prefix + p + suffix for p in pattern_options.split(",")]


def filter_index(entries: dict, pattern: str = "*"):
    """Filter entries by the glob-expanded pattern ``pattern``."""
    if pattern == "*":
        return entries
    patterns = _expand_pattern(pattern)

    ret = {}
    for key, value in entries.items():
        if any(fnmatch.fnmatch(key, pattern) for pattern in patterns):
            ret[key] = value
    return ret


def sort_index(entries: dict, sort_key: str):
    """Sort a dictionary by a existing key in the entries of the dictionary."""
    try:
        return dict(sorted(entries.items(), key=lambda v: getattr(v[1], sort_key)))
    except AttributeError as e:
        raise AttributeError(
            f"The sort key {sort_key} is not a valid attribute of the NNPDF observables."
        )
