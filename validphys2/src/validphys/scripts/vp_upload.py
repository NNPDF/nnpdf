"""
Upload a resource (but not a fit) to the NNPDF server. By default this uploads
a validphys report to the server. For that use `vp-upload <output folder>`.
NB: to upload a fit, use the separate script `vp-uploadfit`.
"""
#Note that the imports are done as late as possible to improve the speed of
#the command line.

import sys


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', help="Folder to upload.")
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


    from validphys import uploadutils
    uploader = uploadutils.ReportUploader()
    try:
        with uploader.upload_or_exit_context(output):
            pass
    except KeyboardInterrupt:
        print(colors.t.bold_red("\nInterrupted by user. Exiting."), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()
