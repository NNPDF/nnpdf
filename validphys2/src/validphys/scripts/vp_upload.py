"""
Simple utility to upload some folder to the NNPDF server.
"""
#Note that the imports are done as late as possible to improve the speed of
#the command line.

import sys


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Upload output to the NNPDF server.")
    parser.add_argument("output", help="Folder to upload.")
    args = parser.parse_args()
    output = args.output

    import os.path as osp
    import logging
    from reportengine import colors
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.addHandler(colors.ColorHandler())

    if not osp.isdir(output):
        log.error("Not a directory: %s", output)
        sys.exit(1)


    from validphys.uploadutils import upload_or_exit_context
    try:
        with upload_or_exit_context(output):
            pass
    except KeyboardInterrupt:
        print(colors.t.bold_red("\nInterrupted by user. Exiting."), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()
