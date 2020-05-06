"""
vp-list

Script which lists available resources locally and remotely

"""
import argparse
import sys
import inspect
import re
import logging

from reportengine import colors

from validphys.loader import FallbackLoader as L
from validphys.config import ConfigError

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())


def atoi(text):
    """convert string to integer, if possible"""
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    sort a list according to natural ordering
    http://nedbatchelder.com/blog/200712/human_sorting.html

    taken directly from https://stackoverflow.com/a/5967539

    """
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "resource",
        type=str,
        help=(
            "The type of resource to check availability for "
            "(locally and remotely). Alternatively type `list` to see what "
            "resources can be checked"
        ),
    )
    args = parser.parse_args()
    l = L()
    attrs = inspect.getmembers(l)
    available = {
        key.lstrip("available_"): val
        for (key, val) in attrs
        if key.startswith("available")
    }
    downloadable = {
        key.lstrip("downloadable_"): val
        for (key, val) in attrs
        if key.startswith("downloadable")
    }
    if args.resource == "list":
        log.info(
            "You can check the availability (locally and remotely) of the following resources:\n- "
            + "\n- ".join(list({*available.keys(), *downloadable.keys()}))
        )
        sys.exit(0)

    if args.resource not in available and args.resource not in downloadable:
        e = ConfigError(
            f"Couldn't find that type of resource in available or downloadable resources.",
            bad_item=args.resource,
            alternatives={*available.keys(), *downloadable.keys()},
        )
        log.error(e)
        print("Alternatively run `vp-list list` to see resources which can be checked.")
        sys.exit(1)
    if args.resource in available and available[args.resource]:
        log.info(
            f"The following {args.resource} are available locally:\n- "
            + "\n- ".join(sorted(available[args.resource], key=natural_keys))
        )
    if args.resource in downloadable and downloadable[args.resource]:
        log.info(
            f"The following {args.resource} are downloadable:\n- "
            + "\n- ".join(sorted(downloadable[args.resource], key=natural_keys))
        )


if __name__ == "__main__":
    main()
