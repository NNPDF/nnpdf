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
import argparse
import pathlib
import logging
import hashlib

from reportengine import app
from validphys.config import Environment, EnvironmentError_
from validphys.app import App

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
        super(SetupFitEnvironment, self).__init__(output='', **kwargs)

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


class SetupFitApp(App):
    """The class which parsers and perform the filtering"""
    environment_class = SetupFitEnvironment

    def __init__(self):
        super(SetupFitApp, self).__init__(name='setup-fit',
                                          providers=['validphys.plots'])

    @property
    def default_style(self):
        return str(self.this_folder() / '../small.mplstyle')

    @property
    def argparser(self):
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('runcard', help="configuration file name")

        cout = parser.add_mutually_exclusive_group()
        cout.add_argument('--cout', action='store_const', const=True,
                          help="display C output. Default depends on log level")
        cout.add_argument('--no-cout', dest='cout',
                          action='store_const', const=False)

        loglevel = parser.add_mutually_exclusive_group()
        loglevel.add_argument('-q', '--quiet', help="supress INFO messages",
                              action='store_true')
        loglevel.add_argument('-d', '--debug', help="show debug info",
                              action='store_true')

        parser.add_argument('-x', '--extra-providers', nargs='+',
                            help="additional providers from which to "
                                 "load actions. Must be an importable specifiaction.")

        parallel = parser.add_mutually_exclusive_group()
        parallel.add_argument('--parallel', action='store_true',
                              help="execute actions in parallel")
        parallel.add_argument('--no-parrallel', dest='parallel',
                              action='store_false')
        return parser

    def run(self):
        try:
            self.environment.config_yml = pathlib.Path(self.args['runcard']).absolute()
            self.environment.init_output()

            # very ugly, require changes/hooks in validphys.app.run
            self.args['config_yml'] = self.environment.config_yml
            self.args['upload'] = False
            self.args['output'] = ''

            # proceed with default run
            super().run()

            # if succeeded print md5
            self.environment.save_md5()
        except SetupFitError as e:
            log.error(f"Error in setup-fit:\n{e}")
            sys.exit(1)
        except Exception as e:
            log.critical(f"Bug in setup-fit ocurred. Please report it.")
            raise

    def get_config(self):
        # first read the runcard
        config_file = self.args['config_yml']
        #return self.config_class.from_yaml(, environment=self.environment)



def main():
    a = SetupFitApp()
    a.main()
