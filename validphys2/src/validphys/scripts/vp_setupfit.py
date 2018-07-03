#!/usr/bin/env python
"""
    setup-fit - prepare and apply data cuts before fit
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
        if not extension:
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

    @classmethod
    def ns_dump_description(cls):
        return {'filter_path': "The filter output folder",
                **super().ns_dump_description()}


class SetupFitConfig(Config):
    """Specialization for yaml parsing"""
    @classmethod
    def from_yaml(cls, o, *args, **kwargs):
        try:
            file_content = yaml.safe_load(o)
        except yaml.error.YAMLError as e:
            raise ConfigError(f"Failed to parse yaml file: {e}")
        file_content['theoryid'] = {'from_': 'theory'}
        file_content['use_cuts'] = False
        file_content['t0pdfset'] = {'from_': 'datacuts'}
        file_content['combocuts'] = {'from_': 'datacuts'}
        file_content['q2min'] = {'from_': 'datacuts'}
        file_content['w2min'] = {'from_': 'datacuts'}
        file_content['rngalgo'] = {'from_': 'fitting'}
        file_content['seed'] = {'from_': 'fitting'}
        file_content['fakedata'] = {'from_': 'closuretest'}
        file_content['fakenoise'] = {'from_': 'closuretest'}
        file_content['fakepdf'] = {'from_': 'closuretest'}
        file_content['filterseed'] = {'from_': 'closuretest'}
        file_content['rancutmethod'] = {'from_': 'closuretest'}
        file_content['rancutprob'] = {'from_': 'closuretest'}
        file_content['errorsize'] = {'from_': 'closuretest'}
        file_content['rancuttrnval'] = {'from_': 'closuretest'}
        file_content['actions_'] = ['check_positivity', 'filter']
        return cls(file_content, *args, ** kwargs)


class SetupFitApp(App):
    """The class which parsers and perform the filtering"""
    environment_class = SetupFitEnvironment
    config_class = SetupFitConfig

    def __init__(self):
        super(SetupFitApp, self).__init__(name='setup-fit',
                                          providers=['validphys.filters'])

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
