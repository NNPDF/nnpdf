"""
uploadutils.py

Tools to upload resources to remote servers.
"""
import subprocess
import logging
import shutil
import uuid
import base64
import sys
import contextlib
import pathlib
import tempfile
from urllib.parse import urljoin

from reportengine.compat import yaml
from reportengine.colors import t

from NNPDF import get_profile_path

log = logging.getLogger(__name__)

class UploadError(Exception): pass

class BadSSH(UploadError): pass

def _profile_key(k):
    """Return a property that fetches a given key from ``self._profile``."""
    @property
    def f(self):
        try:
            return self._profile[k]
        except KeyError as e:
            raise UploadError(f"Profile '{get_profile_path()}' does not contain key '{k}'") from e

    return f

class Uploader():
    """Base class for implementing upload behaviour. The main abstraction is a
    context manager ``upload_context`` which checks that the upload seems
    possible, then does the work inside the context and then uploads the
    result. The various derived classes should be used."""

    def __init__(self):
        self._lazy_profile = None

    upload_host = _profile_key('upload_host')

    @property
    def _profile(self):
        if self._lazy_profile is None:
            with open(get_profile_path()) as f:
                self._lazy_profile = yaml.safe_load(f)
        return self._lazy_profile

    def get_relative_path(self, output_path):
        """Return the relative path to the ``target_dir``."""
        return base64.urlsafe_b64encode(uuid.uuid4().bytes).decode()

    def check_auth(self):
        """Check that we can authenticate with a certificate."""
        ssh_command_line = ('ssh', '-o', 'PreferredAuthentications=publickey',
                            '-q', self.upload_host, 'exit')

        str_line = ' '.join(repr(ele) for ele in ssh_command_line)

        log.info("Checking SSH connection to %s.", self.upload_host)

        try:
            subprocess.run(ssh_command_line, check=True)
        except subprocess.CalledProcessError as e:
            raise BadSSH(("Could not validate the SSH key. "
            "The command\n%s\nreturned a non zero exit status. "
            "Please make sure thet your public SSH key is on the server.")
            % str_line) from e
        except OSError as e:
            raise BadSSH("Could not run the command\n%s\n: %s" % (str_line, e)) from e

        log.info("Connection seems OK.")


    def check_rsync(self):
        """Check that the rsync command exists"""
        if not shutil.which('rsync'):
            raise BadSSH("Could not find the rsync command. "
            "Please make sure it is installed.")


    def upload_output(self, output_path):
        """Rsync ``output_path`` to the server and print the resulting URL. If
        specific_file is given"""
        #Set the date to now
        pathlib.Path(output_path).touch()
        randname = self.get_relative_path(output_path)
        newdir = self.target_dir + randname

        rsync_command = ('rsync', '-aLz', '--chmod=ug=rwx,o=rx',
                         f"{output_path}/", f'{self.upload_host}:{newdir}')

        log.info(f"Uploading output ({output_path}) to {self.upload_host}")
        try:
            subprocess.run(rsync_command, check=True)
        except subprocess.CalledProcessError as e:
            msg = f"Failed to upload output: {e}"
            raise BadSSH(msg) from e
        return randname


    def _print_output(self, name):
        url = urljoin(self.root_url, name)
        log.info(f"Upload completed. The result is available at:\n{t.bold_blue(url)}")


    def check_upload(self):
        """Check that it looks possible to upload something.
        Raise an UploadError if not."""
        self.check_rsync()
        self.check_auth()

    @contextlib.contextmanager
    def upload_context(self, output):
        """Before entering the context, check that uploading is feasible.
        On exiting the context, upload output.
        """
        self.check_upload()
        yield
        res = self.upload_output(output)
        self._print_output(res)

    @contextlib.contextmanager
    def upload_or_exit_context(self, output):
        """Like upload context, but log and sys.exit on error"""
        try:
            with self.upload_context(output):
                yield
        except BadSSH as e:
            log.error(e)
            sys.exit()


class ReportUploader(Uploader):
    """An uploader for validphys reports."""
    target_dir = _profile_key('reports_target_dir')
    root_url = _profile_key('reports_root_url')



class FileUploader(Uploader):
    """Uploader for individual files for single-file resources. It does the "
    "same but prints the URL of the file."""
    def _print_output(self, result, name):
        url = urljoin(result, name)
        log.info(f"Upload completed. The result is available at:\n{t.bold_blue(url)}")

    @contextlib.contextmanager
    def upload_context(self, output_and_file):
        output, specific_file = output_and_file
        self.check_upload()
        yield
        res = self.upload_output(output)
        self._print_output(self.root_url+'/'+res+'/', specific_file)

class ReportFileUploader(FileUploader, ReportUploader):
    pass

class FitUploader(FileUploader):
    """An uploader for fits. Fits will be automatically compressed
    before uploading."""
    target_dir = _profile_key('fits_target_dir')
    root_url = _profile_key('fits_root_url')

    def get_relative_path(self, output_path=None):
        return ''

    def compress(self, output_path):
        """Compress the folder and put in in a directory inside its parent."""
        #make_archive fails if we give it relative paths for some reason
        output_path = output_path.resolve()
        tempdir = tempfile.mkdtemp(prefix='fit_upload_deleteme_',
                                   dir=output_path.parent)
        log.info(f"Compressing fit to {tempdir}")
        archive_path_without_extension = pathlib.Path(tempdir)/(output_path.name)
        try:
            shutil.make_archive(base_name=archive_path_without_extension,
                                format='gztar',
                            root_dir=output_path.parent, base_dir=output_path.name)
        except Exception as e:
            log.error(f"Couldn't compress archive: {e}")
            raise UploadError(e) from e
        return tempdir, archive_path_without_extension


    def upload_output(self, output_path):
        output_path = pathlib.Path(output_path)
        new_out, name = self.compress(output_path)
        super().upload_output(new_out)

        shutil.rmtree(new_out)
        return name.with_suffix('.tar.gz').name

    @contextlib.contextmanager
    def upload_context(self, output_path):
        self.check_upload()
        yield
        res = self.upload_output(output_path)
        self._print_output(self.root_url, res)
