"""
Upload a resource to the NNPDF server.

The script automatically detects (:py:func:`validphys.uploadutils.check_input`)
the type of the input.

 - A ``fit`` is defined to be any folder structure that contains a ``filter.yml`` file at its root
 - A ``hyperscan`` is a ``fit`` that contains ``tries.json`` file without a ``postfit`` folder.
 - a ``PDF`` is any folder containing a ``.info`` file at the root and a replica 0
 - a report is any such structure containing an ``index.html`` file at the root.

The input folder is then placed in the correct location in the server accordingly.

"""
# Note that the imports are done as late as possible to improve the speed of
# the command line.

import sys


def main():
    import argparse

    # Parse the __doc__ str to remove the rtd formatting
    doc_help = __doc__.replace("``", "'")
    doc_help = doc_help.replace("(:py:func:`validphys.uploadutils.check_input`)\n", "")
    parser = argparse.ArgumentParser(
        description=doc_help, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('output', help="Folder to upload.")
    parser.add_argument(
        '-i', '--interactive', help="Interactively create a meta file.", action='store_true'
    )
    parser.add_argument(
        '-f',
        '--force',
        help="If the fit to upload already exists on the server, overwrite it.",
        action='store_true',
    )
    args = parser.parse_args()
    import pathlib

    output = pathlib.Path(args.output)
    interactive = args.interactive
    force = args.force

    import logging
    import os.path as osp

    from reportengine import colors

    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.addHandler(colors.ColorHandler())

    if not osp.isdir(output):
        log.error("Not a directory: %s", output)
        sys.exit(1)

    from validphys import uploadutils

    if interactive:
        uploadutils.interactive_meta(output)

    input_type = uploadutils.check_input(output)
    log.info(f"Detected {input_type} input")

    uploader_dict = {
        'report': uploadutils.ReportUploader,
        'hyperscan': uploadutils.HyperscanUploader,
        'fit': uploadutils.FitUploader,
        'pdf': uploadutils.PDFUploader,
    }
    uploader = uploader_dict[input_type]()

    if isinstance(uploader, uploadutils.ReportUploader):
        upload_args = (output,)
    else:
        upload_args = (output, force)

    try:
        with uploader.upload_or_exit_context(*upload_args):
            pass
    except KeyboardInterrupt:
        print(colors.t.bold_red("\nInterrupted by user. Exiting."), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()
