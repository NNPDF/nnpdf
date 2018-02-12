"""
A more interactive version of vp_upload
"""

#Note that the imports are done as late as possible to improve the speed of
#the command line.

import sys
import pathlib
import os

import logging

import pygments
from prompt_toolkit.shortcuts import prompt

from reportengine import colors


log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

def yes_no_str(default=None):
    """Return a yes or no string for the prompt, with the default
    highlighted"""
    y = 'y' if default is not True else colors.t.bold('Y')
    n = 'n' if default is not False else colors.t.bold('N')
    return f'[{y}/{n}]'



def confirm(message, default=None):
    """
    This is like prompt_toolkit.shortcuts.create_confirm_application
    except that it doesn't bind control+c to "No", which is nonsensical.
    Instead, it raises an exception.

    Also support defaults.
    """
    from prompt_toolkit.enums import DEFAULT_BUFFER
    from prompt_toolkit.key_binding.bindings.basic import load_abort_and_exit_bindings
    from prompt_toolkit.keys import Keys
    from prompt_toolkit.shortcuts import create_prompt_application, run_application

    registry = load_abort_and_exit_bindings()

    @registry.add_binding('y')
    @registry.add_binding('Y')
    def yes_event(event):
        event.cli.buffers[DEFAULT_BUFFER].text = 'y'
        event.cli.set_return_value(True)

    @registry.add_binding('n')
    @registry.add_binding('N')
    def no_event(event):
        event.cli.buffers[DEFAULT_BUFFER].text = 'n'
        event.cli.set_return_value(False)

    #Require explicitly True or False
    if default is True:
        registry.add_binding(Keys.Enter)(yes_event)
    elif default is False:
        registry.add_binding(Keys.Enter)(no_event)
    #There doesn't seem to be an easy way to do this "idiomatically". The lib
    #esccapes ANSI sequences.
    message = f'{message} {yes_no_str(default)}'
    print(message)

    app = create_prompt_application('', key_bindings_registry=registry)
    return run_application(app, patch_stdout=True)

def handle_single_file(filename):
    import tempfile
    out = pathlib.Path(tempfile.mkdtemp(prefix='vp-upload'))
    filename = pathlib.Path(filename)
    p = out / filename.name
    p.symlink_to(filename.absolute())
    return out, filename.name

def edit_settings(d):
    title = d.get('title', '')
    author = d.get('author', '')
    keywords = d.get('keywords', '')

    d['title'] = prompt('title: ', default=title)

    if not author:
        try:
            import pwd
        except ImportError:
            pass
        else:
            author = pwd.getpwuid(os.getuid())[4]

    d['author'] = prompt("author: ", default=author)
    kwinp = prompt("keywords: ", default=','.join(keywords))
    d['keywords'] = [k.strip() for k in kwinp.split(',') if k]

def handle_meta_interactive(output):
    metapath = output / 'meta.yaml'
    import ruamel_yaml as yaml
    #The yaml lexer is broken. Use something else.
    lex = pygments.lexers.get_lexer_by_name('pkgconfig')
    fmt = pygments.formatters.TerminalFormatter()
    if metapath.exists():
        log.info("Found meta.yaml file")

        with open(metapath) as f:
            content = f.read()

        print(pygments.highlight(content, lex, fmt))

        msg = "Use these settings? (answering no will edit the meta.yaml file)"
        edit = not confirm(msg, default=True)

        if edit:
            d = yaml.load(content, yaml.RoundTripLoader)
        else:
            return

    else:
        #We are making these the empty string, because prompt_toolkit doesn't
        #support default=None.
        d = {'title': '', 'author': '', 'keywords':''}

    import io
    while True:
        edit_settings(d)

        print("Metadata:")

        s = io.StringIO()
        yaml.dump(d, s, yaml.RoundTripDumper)
        metastr = s.getvalue()
        print(pygments.highlight(metastr, lex, fmt))

        if confirm("Confirm?"):
            break


    with open(metapath, 'w') as f:
        f.write(metastr)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Upload output to the NNPDF server.")
    parser.add_argument("output", help="Folder to upload.")
    args = parser.parse_args()
    output = pathlib.Path(args.output)
    upload_output = output

    from validphys import uploadutils

    if not output.is_dir():
        if output.is_file():
            upargs = handle_single_file(output)
            upload_output = upargs[0]
            uploader = uploadutils.SingleReportFileUploader()
        else:
            if not output.exists():
                log.error(f"No such file or directory: {output}")
            else:
                log.errorr(f"{output} is not a file or a directory")
            sys.exit(1)
    else:
        uploader = uploadutils.ReportUploader()
        upargs = output


    try:
        with uploader.upload_or_exit_context(upargs):
            handle_meta_interactive(upload_output)
    except (KeyboardInterrupt, EOFError):
        print(colors.t.bold_red("\nInterrupted by user. Exiting."), file=sys.stderr)
        exit(1)


if __name__ == '__main__':
    main()