"""
uploadutils.py

Tools to upload resources to remote servers.
"""
import base64
import contextlib
from glob import glob
import hashlib
import logging
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.parse import urljoin
import uuid

import prompt_toolkit
from prompt_toolkit.completion import WordCompleter

from reportengine.colors import t
from reportengine.compat import yaml
from validphys.loader import Loader, RemoteLoader
from validphys.renametools import Spinner

log = logging.getLogger(__name__)


class UploadError(Exception):
    pass


class BadSSH(UploadError):
    pass


def _profile_key(k):
    """Return a property that fetches a given key from ``self._profile``."""

    @property
    def f(self):
        try:
            return self._profile[k]
        except KeyError as e:
            raise UploadError(f"Profile '{self._profile_path}' does not contain key '{k}'") from e

    return f


class Uploader:
    """Base class for implementing upload behaviour. The main abstraction is a
    context manager ``upload_context`` which checks that the upload seems
    possible, then does the work inside the context and then uploads the
    result. The various derived classes should be used."""

    upload_host = _profile_key('upload_host')
    _profile = Loader().nnprofile

    def get_relative_path(self, output_path):
        """Return the relative path to the ``target_dir``."""
        return base64.urlsafe_b64encode(uuid.uuid4().bytes).decode()

    def check_auth(self):
        """Check that we can authenticate with a certificate."""
        ssh_command_line = (
            'ssh',
            '-o',
            'PreferredAuthentications=publickey',
            '-q',
            self.upload_host,
            'exit',
        )

        str_line = ' '.join(repr(ele) for ele in ssh_command_line)

        log.info("Checking SSH connection to %s.", self.upload_host)

        try:
            subprocess.run(ssh_command_line, check=True)
        except subprocess.CalledProcessError as e:
            raise BadSSH(
                (
                    "Could not validate the SSH key. "
                    "The command\n%s\nreturned a non zero exit status. "
                    "Please make sure that your public SSH key is on the server."
                )
                % str_line
            ) from e
        except OSError as e:
            raise BadSSH("Could not run the command\n%s\n: %s" % (str_line, e)) from e

        log.info("Connection seems OK.")

    def check_rsync(self):
        """Check that the rsync command exists"""
        if not shutil.which('rsync'):
            raise BadSSH("Could not find the rsync command. Please make sure it is installed.")

    def upload_output(self, output_path):
        """Rsync ``output_path`` to the server and print the resulting URL. If
        specific_file is given"""
        # Set the date to now
        pathlib.Path(output_path).touch()
        randname = self.get_relative_path(output_path)
        newdir = self.target_dir + randname

        rsync_command = (
            'rsync',
            '-aLz',
            '--chmod=ug=rwx,o=rx',
            f"{output_path}/",
            f"{self.upload_host}:{newdir}",
        )

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
        self._print_output(self.root_url + '/' + res + '/', specific_file)


class ReportFileUploader(FileUploader, ReportUploader):
    pass


class ArchiveUploader(FileUploader):
    """Uploader for objects comprising many files such as fits or PDFs"""

    target_dir = None
    root_url = None
    _loader_name = None  # vp loader for this kind of archive
    _resource_type = "Archive"  # name used during logging

    def get_relative_path(self, output_path=None):
        return ''

    def _check_existence(self, resource_name):
        """Check whether the given resource exists on the server.
        Returns true if the resource exists with the same name on the server
        or false otherwise.
        Note that the type of resource being checked is defined by the ``_loader_name`` attribute
        """
        l = RemoteLoader()
        resource_list = getattr(l, self._loader_name)

        return resource_name in resource_list

    def _check_is_indexed(self, resource_name):
        """Check whether the fit is correctly indexed in the server"""
        log.info("Checking whether %s was correctly uploaded...", resource_name)
        time.sleep(3)
        if self._check_existence(resource_name):
            log.info("It has been correctly indexed by the server!")
        else:
            log.error(
                "The object is uploaded but hasn't been indexed yet by the server. "
                "You should upload it again to ensure it is indexed: vp-upload %s",
                resource_name,
            )

    def _compress(self, output_path):
        """Compress the folder and put in in a directory inside its parent."""
        # make_archive fails if we give it relative paths for some reason
        output_path = output_path.resolve()
        tempdir = tempfile.mkdtemp(
            prefix=f'{self._resource_type}_upload_deleteme_', dir=output_path.parent
        )
        log.info(f"Compressing {self._resource_type} to {tempdir}")
        archive_path_without_extension = pathlib.Path(tempdir) / (output_path.name)
        try:
            with Spinner():
                shutil.make_archive(
                    base_name=archive_path_without_extension,
                    format='gztar',
                    root_dir=output_path.parent,
                    base_dir=output_path.name,
                )
        except Exception as e:
            log.error(f"Couldn't compress archive: {e}")
            raise UploadError(e) from e
        return tempdir, archive_path_without_extension

    def upload_output(self, output_path, force):
        output_path = pathlib.Path(output_path)
        fit_name = output_path.name

        if not force:
            if self._check_existence(fit_name):
                log.error(
                    "A %s with the same name already exists on "
                    "the server. To overwrite it use the "
                    "--force flag, as in `vp-upload <%s_name> --force.",
                    self._resource_type,
                    self._resource_type,
                )
                raise UploadError

        new_out, name = self._compress(output_path)
        super().upload_output(new_out)

        # Check whether the fit was really uploaded
        self._check_is_indexed(fit_name)

        shutil.rmtree(new_out)
        return name.with_suffix('.tar.gz').name

    @contextlib.contextmanager
    def upload_context(self, output_path, force):
        self.check_upload()
        yield
        res = self.upload_output(output_path, force)
        self._print_output(self.root_url, res)

    @contextlib.contextmanager
    def upload_or_exit_context(self, output, force):
        try:
            with self.upload_context(output, force):
                yield
        except BadSSH as e:
            log.error(e)
            sys.exit()


class FitUploader(ArchiveUploader):
    """An uploader for fits. Fits will be automatically compressed
    before uploading."""

    target_dir = _profile_key('fits_target_dir')
    root_url = _profile_key('fits_root_url')
    _loader_name = "downloadable_fits"
    _resource_type = "fit"

    def check_fit_md5(self, output_path):
        """When ``vp-setupfit`` is run successfully, it creates an ``md5`` from
        the config. We check that the ``md5`` matches the ``filter.yml`` which
        is checking that ``vp-setupfit`` was run and that the ``filter.yml``
        inside the fit folder wasn't modified.

        """
        md5_path = output_path / "md5"
        try:
            with open(md5_path, "r") as f:
                saved_md5 = f.read()
        except FileNotFoundError as e:
            log.error(
                "It doesn't appear that `vp-setupfit` was run because no `md5` "
                "was found, `vp-setupfit` should be run before uploading a fit."
            )
            raise UploadError(f"Fit MD5 file not found at {md5_path}") from e

        with open(output_path / "filter.yml", "rb") as f:
            hashed_config = hashlib.md5(f.read()).hexdigest()

        if hashed_config != saved_md5:
            log.error(
                "Saved md5 doesn't match saved fit configuration runcard, which "
                "suggests that the configuration file was modified after it was "
                "saved. <fit folder>/filter.yml shouldn't be modified directly. "
                "Instead modify the fit runcard and re-run ``vp-setupfit``."
            )
            raise UploadError

    def upload_output(self, output_path, force):
        output_path = pathlib.Path(output_path)
        self.check_fit_md5(output_path)
        return super().upload_output(output_path, force)


class HyperscanUploader(FitUploader):
    """Uploader for hyperopt scans, which are just special cases of fits"""

    _resource_type = "hyperscans"
    _loader_name = "downloadable_hyperscans"
    target_dir = _profile_key('hyperscan_target_dir')
    root_url = _profile_key('hyperscan_root_url')


class PDFUploader(ArchiveUploader):
    """An uploader for PDFs. PDFs will be automatically compressed
    before uploading."""

    target_dir = _profile_key('pdfs_target_dir')
    root_url = _profile_key('pdfs_root_url')
    _loader_name = "downloadable_pdfs"
    _resource_type = "pdf"


def check_for_meta(path):
    """Function that checks if a report input has a ``meta.yaml`` file.
    If not it prompts the user to either create one or follow an interactive
    prompt which assists the user in creating one.

    Parameters
    ----------
    path: pathlib.Path
        Input path

    Returns
    -------
    None
    """
    if "meta.yaml" not in os.listdir(path):
        raise FileNotFoundError(
            "No meta.yaml file found. Please either add "
            "the meta tags to the runcard or use the --interactive flag "
            "with vp-upload to interactively create one"
        )
    return True


def interactive_meta(path):
    """Function to interactively create a meta.yaml file

    Parameters
    ----------
    path: pathlib.Path
        Input path

    Returns
    -------
    None
    """
    # Import here to avoid circular imports
    from validphys.promptutils import KeywordsWithCache

    title = prompt_toolkit.prompt("Enter report title: ")

    default = ""
    try:
        import pwd
    except ImportError:
        pass
    else:
        default = pwd.getpwuid(os.getuid())[0]
    author = prompt_toolkit.prompt("Enter author name: ", default=default)

    kwinp = prompt_toolkit.prompt(
        "Enter keywords: ",
        completer=WordCompleter(words=KeywordsWithCache(RemoteLoader())),
        complete_in_thread=True,
    )
    keywords = [k.strip() for k in kwinp.split(",") if k]

    meta_dict = {"title": title, "author": author, "keywords": keywords}
    with open(path / "meta.yaml", "w") as stream:
        yaml.safe_dump(meta_dict, stream)


def check_input(path):
    """A function that checks the type of the input for vp-upload. The type
    determines where on the vp server the file will end up

    A ``fit`` is defined as any folder structure containing a ``filter.yml``
    file at its root.

    A ``pdf`` is defined as any folder structure that contains a ``.info``
    file and a replica 0 at its root.

    A ``report`` is defined as any folder structure that contains an ``index.html``
    at its root.

    If the input file does not fall under any such category ``ValueError`` exception
    is raised and the user is prompted to use either ``rsync`` or
    :py:mod:`validphys.scripts.wiki_upload`.

    Parameters
    ----------
    path: pathlib.Path
        Path of the input file
    """
    files = os.listdir(path)
    # Require that a .info file and replica 0 exist before admitting
    # the input is a valid LHAPDF set
    info_reg, rep0_reg = map(re.compile, ('.+\.info', '.+0000\.dat'))

    if 'meta.yaml' in files:
        return 'report'
    elif 'index.html' in files and check_for_meta(path):
        # It could be that the user has created a report but hasn't
        # created a meta.yaml file. In which case we raise an exception
        # and instruct the user to either create one or use the
        # --interactive flag to create one
        return 'report'
    elif 'filter.yml' in files:
        # The product of a n3fit run, usually a fit but could be a hyperopt scan
        # For that there should be a) tries.json files and b) no postfit
        if "postfit" not in files and glob(path.as_posix() + "/nnfit/replica_*/tries.json"):
            return 'hyperscan'
        return 'fit'
    elif list(filter(info_reg.match, files)) and list(filter(rep0_reg.match, files)):
        return 'pdf'
    else:
        log.error(
            f"Specified input directory: {path} did not fall under the known "
            "categories of validphys (report, fit, or pdf)."
        )
        raise ValueError(
            "Unrecognized type of input, "
            "please save to the server using rsync or wiki-upload. "
            "The --interactive flag will generate a meta file which "
            "will cause the input to be registered as a report."
        )
