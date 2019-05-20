#!/usr/bin/env python
"""
    n3fit - performs fit using ml external frameworks
"""

import sys
import re
import shutil
import pathlib
import logging
import warnings
import argparse

from validphys.app import App
from validphys.config import Environment, Config
from validphys.config import EnvironmentError_, ConfigError
from reportengine import colors
from reportengine.compat import yaml


N3FIT_FIXED_CONFIG = dict(
    use_cuts = 'internal',
    use_t0 = True,
    actions_ = ['datacuts::theory fit']
)

N3FIT_PROVIDERS = ['fit']

log = logging.getLogger(__name__)

RUNCARD_COPY_FILENAME = 'filter.yml'


class N3FitError(Exception):
    """Exception raised when n3fit cannot succeed and knows why"""
    pass


class N3FitEnvironment(Environment):
    """Container for information to be filled at run time"""
    def init_output(self):
        # check file exists, is a file, has extension.
        if not self.config_yml.exists():
            raise N3FitError("Invalid runcard. File not found.")
        else:
            if not self.config_yml.is_file():
                raise N3FitError("Invalid runcard. Must be a file.")

        # check if results folder exists
        self.output_path = pathlib.Path(self.output_path).absolute()
        if not self.output_path.exists():
            if not re.fullmatch(r'[\w.\-]+', self.output_path.name):
                raise N3FitError("Invalid output folder name. Must be alphanumeric.")
            try:
                self.output_path.mkdir()
                pathlib.Path(self.output_path / 'nnfit').mkdir()
            except OSError as e:
                raise EnvironmentError_(e) from e

            try:
                shutil.copy2(self.config_yml, self.output_path / RUNCARD_COPY_FILENAME)
            except shutil.SameFileError:
                pass

        # create output folder for the fit
        self.replica_path = self.output_path / 'nnfit'
        for replica in self.replica:
            path = self.replica_path / 'replica_{0}'.format(replica)
            log.info('Creating replica output folder in {0}'.format(path))
            try:
                path.mkdir(exist_ok=True)
            except OSError as e:
                raise EnvironmentError_(e) from e


    @classmethod
    def ns_dump_description(cls):
        return {'replica': 'The MC replica number',
                'replica_path': 'The replica output path',
                'output_path': 'The runcard name',
                'hyperopt': 'The hyperopt flag',
                'create_test_card' : 'Create test card flag',
                **super().ns_dump_description()}


class N3FitConfig(Config):
    """Specialization for yaml parsing"""

    @classmethod
    def from_yaml(cls, o, *args, **kwargs):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',
                                      yaml.error.MantissaNoDotYAML1_1Warning)
                #We need to specify the older version 1.1 to support the
                #older configuration files, which liked to use on/off for
                #booleans.
                #The floating point parsing yields warnings everywhere, which
                #we suppress.
                file_content = yaml.safe_load(o, version='1.1')
        except yaml.error.YAMLError as e:
            raise ConfigError(f"Failed to parse yaml file: {e}")
        if not isinstance(file_content, dict):
            raise ConfigError(f"Expecting input runcard to be a mapping, "
                              f"not '{type(file_content)}'.")
        file_content.update(N3FIT_FIXED_CONFIG)
        return cls(file_content, *args, **kwargs)


class N3FitApp(App):
    """The class which parsers and performs the fit"""
    environment_class = N3FitEnvironment
    config_class = N3FitConfig

    def __init__(self):
        super(N3FitApp, self).__init__(name='n3fit',
                                       providers=N3FIT_PROVIDERS)

    @property
    def argparser(self):
        parser = super().argparser
        parser.add_argument('-o', '--output',
                            help='Output folder and name of the fit',
                            default=None)
        parser.add_argument('--hyperopt',
                            help='Enable hyperopt scan',
                            default=None, type=int)
        parser.add_argument('--create_test_card',
                            help='Generates a new runcard where some datasets are only used for a test set (only hyperopt)',
                            action='store_true')

        def check_positive(value):
            ivalue = int(value)
            if ivalue <= 0:
                raise argparse.ArgumentTypeError("%s is an invalid positive int value." % value)
            return ivalue

        parser.add_argument('replica',
                            help='MC replica number', type=check_positive)
        parser.add_argument('-r', '--replica_range',
                help = 'End of the range of replicas to compute', type=check_positive)
        return parser

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        if args['output'] is None:
            args['output'] = pathlib.Path(args['config_yml']).stem
        return args

    def run(self):
        try:
            self.environment.config_yml = pathlib.Path(self.args['config_yml']).absolute()
            replica = self.args['replica']
            if self.args['replica_range']:
                replicas = list(range(replica, self.args['replica_range']+1))
            else:
                replicas = [replica]
            self.environment.replica = replicas
            self.environment.hyperopt = self.args['hyperopt']
            if self.args['create_test_card']:
                self.environment.create_test_card = self.environment.config_yml
            else:
                self.environment.create_test_card = None
            super().run()
        except N3FitError as e:
            log.error(f'Error in n3fit:\n{e}')
            sys.exit(1)
        except Exception as e:
            log.critical(f'Bug in n3fit ocurred. Please report it.')
            print(
                colors.color_exception(e.__class__, e, e.__traceback__),
                file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    a = N3FitApp()
    a.main()
