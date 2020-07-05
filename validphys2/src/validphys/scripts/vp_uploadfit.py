"""
Upload a fit to the NNPDF server. To do this use `vp-uploadfit <fit folder>`.
"""
#Note that the imports are done as late as possible to improve the speed of
#the command line.

import sys


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', help="Fit folder to upload.")
    parser.add_argument(
        '-f',
        '--force',
        help="If the fit to upload already exists on the server, overwrite it.",
        action='store_true')
    args = parser.parse_args()
    output = args.output
    force = args.force

    import os.path as osp
    import logging
    from reportengine import colors
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.addHandler(colors.ColorHandler())

    if not osp.isdir(output):
        log.error("Not a directory: %s", output)
        sys.exit(1)


    from validphys import uploadutils
    uploader = uploadutils.FitUploader()
    try:
        with uploader.upload_or_exit_context(output, force):
            pass
    except KeyboardInterrupt:
        print(colors.t.bold_red("\nInterrupted by user. Exiting."), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()
