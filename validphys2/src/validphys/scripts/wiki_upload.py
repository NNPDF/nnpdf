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
from prompt_toolkit import HTML
from prompt_toolkit.shortcuts import prompt

from reportengine import colors

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

def yes_no_str(default=None):
    """Return a yes or no string for the prompt, with the default
    highlighted"""
    if default is None:
        return f'[y/n]'
    elif default:
        return HTML('[<b>Y</b>/n]')
    else:
        return HTML('[y/<b>N</b>]')


def confirm(message, default=None):
    """
    This is like prompt_toolkit.shortcuts.confirm (implemented by
    create_confirm_session) except that it doesn't bind control+c to "No", but
    instead raises an exception.

    It also support defaults.
    """
    from prompt_toolkit.key_binding.key_bindings import KeyBindings
    from prompt_toolkit.keys import Keys
    from prompt_toolkit.formatted_text import merge_formatted_text
    from prompt_toolkit.shortcuts import PromptSession
    bindings = KeyBindings()

    @bindings.add('y')
    @bindings.add('Y')
    def yes(event):
        session.default_buffer.text = 'y'
        event.app.exit(result=True)

    @bindings.add('n')
    @bindings.add('N')
    def no(event):
        session.default_buffer.text = 'n'
        event.app.exit(result=False)

    @bindings.add(Keys.Any)
    def nothing(event):
        " Disallow inserting other text. "
        pass

    if default:
        bindings.add(Keys.Enter)(yes)
    elif default is not None:
        bindings.add(Keys.Enter)(no)
    else:
        bindings.add(Keys.Enter)(nothing)

    complete_message = merge_formatted_text([message, yes_no_str(default)])
    session = PromptSession(complete_message, key_bindings=bindings)
    return session.prompt()

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
    from reportengine.compat import yaml
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
            uploader = uploadutils.ReportFileUploader()
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
