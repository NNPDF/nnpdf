#!/usr/bin/env python
"""
    setup-fit - prepare and apply data cuts before fit (filter replacement)

    setup-fit constructs the fit [results] folder where data used by nnfit
    will be stored.
"""
import os
import sys
import re
import shutil
import pathlib
import logging
import hashlib

from validphys.config import Environment, Config, EnvironmentError_, ConfigError
from validphys.app import App
from reportengine.compat import yaml
from reportengine import colors

log = logging.getLogger(__name__)

RUNCARD_COPY_FILENAME = "filter.yml"
FILTER_OUTPUT_FOLDER = "filter"
MD5_FILENAME = "md5"


class SetupFitError(Exception):
    """Exception raised when setup-fit cannot succeed and knows why"""
    pass


class SetupFitEnvironment(Environment):
    """Container for information to be filled at run time"""
    def __init__(self, **kwargs):
        super(SetupFitEnvironment, self).__init__(**kwargs)

    def init_output(self):
        # check file exists, is a file, has extension.
        if not self.config_yml.exists():
            raise SetupFitError("Invalid runcard. File not found.")
        else:
            if not re.fullmatch(r'[\w.\-]+', self.config_yml.name):
                raise SetupFitError("Invalid runcard. Must be alphanumeric.")
            if not self.config_yml.is_file():
                raise SetupFitError("Invalid runcard. Must be a file.")

        # check filename
        filename, extension = os.path.splitext(self.config_yml.name)
        if not len(extension):
            raise SetupFitError("Invalid runcard. File extension missing.")

        # check if results folder exists
        self.output_path = pathlib.Path(filename).absolute()
        if self.output_path.exists():
            log.warning(f"Output folder exists: {self.output_path} Overwritting contents")
        else:
            try:
                self.output_path.mkdir()
            except OSError as e:
                raise EnvironmentError_(e) from e

        try:
            shutil.copy2(self.config_yml, self.output_path / RUNCARD_COPY_FILENAME)
        except shutil.SameFileError:
            pass

        # create output folder
        self.filter_path = self.output_path / FILTER_OUTPUT_FOLDER
        self.filter_path.mkdir(exist_ok=True)

    def save_md5(self):
        """Generate md5 key from file"""
        output_filename = self.output_path / MD5_FILENAME
        with open(self.config_yml, 'rb') as f:
            hash_md5 = hashlib.md5(f.read()).hexdigest()
        with open(output_filename, 'w') as g:
            g.write(hash_md5)
        log.info(f"md5 {hash_md5} stored in {output_filename}")


class SetupFitConfig(Config):
    """Specialization for yaml parsing"""
    @classmethod
    def from_yaml(cls, o, *args, **kwargs):
        try:
            file_content = yaml.round_trip_load(o)
            file_content['use_cuts'] = False
            file_content['actions_'] = ['report(main=true)']
        except yaml.error.YAMLError as e:
            raise ConfigError(f"Failed to parse yaml file: {e}")
        return cls(file_content, *args, ** kwargs)


class SetupFitApp(App):
    """The class which parsers and perform the filtering"""
    environment_class = SetupFitEnvironment
    config_class = SetupFitConfig

    def __init__(self):
        super(SetupFitApp, self).__init__(name='setup-fit',
                                          providers=['validphys.plots'])

    @property
    def default_style(self):
        return str(self.this_folder() / '../small.mplstyle')

    def run(self):
        try:
            # set folder output name
            self.environment.config_yml = pathlib.Path(self.args['config_yml']).absolute()

            # proceed with default run
            super().run()

            # if succeeded print md5
            self.environment.save_md5()
        except SetupFitError as e:
            log.error(f"Error in setup-fit:\n{e}")
            sys.exit(1)
        except Exception as e:
            log.critical(f"Bug in setup-fit ocurred. Please report it.")
            print(
                colors.color_exception(e.__class__, e, e.__traceback__),
                file=sys.stderr)
            sys.exit(1)


def main():
    a = SetupFitApp()
    a.main()
