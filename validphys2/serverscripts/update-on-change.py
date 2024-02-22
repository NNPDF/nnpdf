#!/usr/bin/env python

import argparse

from asyncwatch import EVENTS, watch
import curio
from curio import subprocess


async def main(folder, exe_args):
    watch_events = EVENTS.CREATE | EVENTS.DELETE | EVENTS.MODIFY
    while True:
        async with watch(folder, watch_events):
            await subprocess.run(exe_args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder')
    parser.add_argument('exe_args', nargs='+')
    args = parser.parse_args()
    curio.run(main(args.folder, args.exe_args))
