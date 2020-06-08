"""
vp-list

Script which lists available resources locally and remotely

"""
import argparse
from functools import partial
import re
import logging

from reportengine import colors

from validphys.loader import FallbackLoader as L

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

REMOTE_TOKEN = "downloadable_"
LOCAL_TOKEN = "available_"


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


sane_order = partial(sorted, key=natural_keys)


def main(command_line=None):
    parser = argparse.ArgumentParser(description=__doc__)

    attrs = dir(L)

    available = [
        attr.lstrip(LOCAL_TOKEN) for attr in attrs if attr.startswith(LOCAL_TOKEN)
    ]
    downloadable = [
        attr.lstrip(REMOTE_TOKEN) for attr in attrs if attr.startswith(REMOTE_TOKEN)
    ]
    # set metavar and print choices in help string - otherwise looks ugly.
    parser.add_argument(
        "resource",
        type=str,
        choices={*available, *downloadable},
        help=(
            "The type of resource to check availability for (locally and/or remotely). "
            + "Choose from: "
            + ", ".join(list({*available, *downloadable}))
            + "."
        ),
        metavar="resource",
    )
    g = parser.add_mutually_exclusive_group()
    g.add_argument(
        "-r",
        "--remote-only",
        dest="remote",
        action="store_true",
        default=False,
        help="Only list remote resources",
    )
    g.add_argument(
        "-l",
        "--local-only",
        dest="local",
        action="store_true",
        default=False,
        help="Only list local resources",
    )
    args = parser.parse_args(command_line)

    l = L()
    if not args.remote:
        local_res = getattr(l, LOCAL_TOKEN + args.resource, None)
        if args.resource in available and local_res:
            log.info("The following %s are available locally:", args.resource)
            print("- " + "\n- ".join(sane_order(local_res)))
        else:
            log.info("No %s are available locally.", args.resource)
    if not args.local:
        remote_res = getattr(l, REMOTE_TOKEN + args.resource, None)
        if args.resource in downloadable and remote_res:
            log.info("The following %s are downloadable:", args.resource)
            print("- " + "\n- ".join(sane_order(remote_res)))
        else:
            log.info("No %s are available to download.", args.resource)


if __name__ == "__main__":
    main()
