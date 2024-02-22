#!/usr/bin/env python3
"""
NNPDF resource downloader. The basic syntax is

vp-get <resource_type> <resource_name>

Use

vp-get --list

to see a list of resource types.
If the resource is already installed, a string with its name will be
printed to stdout. If not, it will be searched in the remote repositories
and installed if found.
"""
import argparse
import logging
import sys

from reportengine import colors
from reportengine.baseexceptions import ErrorWithAlternatives
from validphys.loader import FallbackLoader as Loader
from validphys.loader import LoadFailedError

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())


class ListAction(argparse.Action):
    def __init__(self, *args, loader, **kwargs):
        self.loader = loader
        super().__init__(*args, **kwargs)

    def __call__(self, parser, *args, **kwargs):
        prefix = 'download_'
        prefix_len = len(prefix)
        tps = [f'\n - {it[prefix_len:]}' for it in dir(self.loader) if it.startswith(prefix)]
        print(f"Available resource types:{''.join(tps)}")
        parser.exit()


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    try:
        l = Loader()
    except Exception as e:
        log.error(f"Could not set up environment: {e}")
        sys.exit(1)
    p.add_argument(
        'resource_type',
        help="Type of the resource to be obtained. " "See --list for a list of resource types.",
    )
    p.add_argument('resource_name', help="Identifier of the resource.")
    p.add_argument(
        '--list', action=ListAction, loader=l, nargs=0, help="List available resources and exit."
    )
    args = p.parse_args()

    tp = args.resource_type
    name = args.resource_name

    try:
        f = getattr(l, f'check_{tp}')
    except AttributeError as e:
        sys.exit(f"No such resource {tp}")

    try:
        res = f(name)
    except LoadFailedError as e:
        print(ErrorWithAlternatives(f"Could not find resource ({tp}): '{name}'.", name))
        sys.exit("Failed to download resource.")
    print(repr(res))
