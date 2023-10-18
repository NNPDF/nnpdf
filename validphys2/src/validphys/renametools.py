"""
    A collection of utility functions to handle logistics of
    LHAPDFs and fits. For use by vp-scripts.
"""
import os
import sys
import threading
import time


class Spinner:
    """Context manager to provide a spinning cursor
    while validphys performs some other task silently.

    When exececuted in a TTY, it shows a spinning cursor for the duration of
    the context manager. In non interactive prompts, it prints to stdout at the
    beginning and end.

    Example
    -------
    >>> from validphys.renametools import Spinner
    >>> with Spinner():
    ...     import time
    ...     time.sleep(5)

    """

    def __init__(self, delay=0.1):
        self.spinner_generator = self.spinning_cursor()
        self.delay = delay

    @property
    def interactive(self):
        return os.isatty(sys.stdout.fileno())

    def spinner_task(self):
        while not self.event.isSet():
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()

    def __enter__(self):
        if self.interactive:
            self.event = threading.Event()
            threading.Thread(target=self.spinner_task).start()
        else:
            print("Waiting...")

    def __exit__(self, exception, value, tb):
        if self.interactive:
            self.event.set()
            sys.stdout.write('\r')
            sys.stdout.flush()
        else:
            print("Done")

    @staticmethod
    def spinning_cursor():
        while True:
            for cursor in '|/-\\':
                yield cursor


def rename_pdf(pdf_folder, initial_fit_name, final_name):
    for item in os.listdir(pdf_folder):
        p = pdf_folder / item
        if p.is_symlink():
            replica = p.resolve().parent.name
            pointer = f'../../nnfit/{replica}/{final_name}.dat'
            p.unlink()
            p.symlink_to(pointer)
        newname = p.name.replace(initial_fit_name, final_name)
        p.rename(p.with_name(newname))
    pdf_folder.rename(pdf_folder.with_name(final_name))


def rename_nnfit(nnfit_path, initial_fit_name, final_name):
    info_file = nnfit_path / f'{initial_fit_name}.info'
    info_file.rename(info_file.with_name(f'{final_name}.info'))
    # Some older fits have the PDF here
    pdf_folder = nnfit_path / initial_fit_name
    if pdf_folder.is_dir():
        rename_pdf(pdf_folder, initial_fit_name, final_name)
    # Change replica names
    for item in nnfit_path.glob('replica*'):
        if item.is_dir():
            files = item.glob(initial_fit_name + '*')
            for i in files:
                newname = i.name.replace(initial_fit_name, final_name)
                i.rename(item / newname)


def rename_postfit(postfit_path, initial_fit_name, final_name):
    pdf_folder = postfit_path / initial_fit_name
    rename_pdf(pdf_folder, initial_fit_name, final_name)
    os.system(f'sed -i -e "s/{initial_fit_name}/{final_name}/g" {postfit_path/"postfit.log"}')


def change_name(initial_path, final_name):
    """Function that takes initial fit name and final fit name
    and performs the renaming"""
    initial_fit_name = initial_path.name
    nnfit = initial_path / 'nnfit'
    if nnfit.exists():
        rename_nnfit(nnfit, initial_fit_name, final_name)
    postfit = initial_path / 'postfit'
    if postfit.exists():
        rename_postfit(postfit, initial_fit_name, final_name)
    newpath = initial_path.with_name(final_name)
    initial_path.rename(newpath)
    return newpath
