# -*- coding: utf-8 -*-
"""
Simple utility to upload some folder to the NNPDF server.

@author: Zahari Kassabov
"""
import sys
import logging

from validphys.app import App

log = logging.getLogger(__name__)

def main():
    import argparse
    import os.path as osp
    parser = argparse.ArgumentParser(description="Upload output to the NNPDF server.")
    parser.add_argument("output", help="Folder to upload.")

    a = App()
    app_args = {'loglevel':logging.INFO}
    a.init_logging(app_args)

    args = parser.parse_args()
    output = args.output
    if not osp.isdir(output):
        log.error("Not a directory: %s", output)
        sys.exit(1)


    with a.upload_context(True, output):
        pass

if __name__ == '__main__':
    main()
