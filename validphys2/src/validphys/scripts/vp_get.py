#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Resource downloader
"""
import sys
import argparse
import logging

from reportengine.baseexceptions import ErrorWithAlternatives
from reportengine import colors
from validphys.loader import FallbackLoader as Loader, LoadFailedError



log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('type')
    p.add_argument('name')
    args = p.parse_args()
    tp = args.type
    name = args.name

    l = Loader()
    try:
        f = getattr(l, f'check_{tp}')
    except AttributeError as e:
        sys.exit(f"No such resource {tp}")

    try:
        res = f(name)
    except LoadFailedError as e:
        raise ErrorWithAlternatives(f"Could not find resource ({tp}): '{name}'.", name)
    print(repr(res))



