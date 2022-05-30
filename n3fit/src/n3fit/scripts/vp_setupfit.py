#!/usr/bin/env python
"""
    setup-fit - prepare and apply data cuts before fit
    setup-fit constructs the fit [results] folder where data used by nnfit
    will be stored.
"""

# Implementation notes
#
# This is a validphys-like app in disguise. It takes an nnfit runcard and adds
# a fixed list of actions and some associated resourced to it so as to make it
# a proper validphys runcard. These config options are defined in the
# SETUPFIT_FIXED_CONFIG mapping below. Similarly, defult options are specified
# in SETUPFIT_DEFAULTS.
#
# Extensions to the setup procedure can be implemented by adding suitable
# actions_ to the mapping (making sure that they are executed in the right
# namespace that pulls all the required resources from the fit runcard),
# together with the additional non variable resources required by said actions
# (such as `use_cuts: "internal"`) in the current code. vp-setupfit also gets
# its own provider modules, so you may need to add the modules of your actions
# to SETUPFIT_PROVIDERS.
#
# The state of the output folder must be such that the nnfit code can be run on
# top.


import sys
import re
import shutil
import pathlib
import logging
import hashlib
import warnings

from validphys.config import Environment, Config, EnvironmentError_, ConfigError
from validphys.app import App
from reportengine.compat import yaml
from reportengine import colors


SETUPFIT_FIXED_CONFIG = dict(
    actions_=[
        'datacuts check_t0pdfset',
        'theory check_positivity',
    ])

SETUPFIT_PROVIDERS = ['validphys.filters',
                      'validphys.theorycovariance.construction',
                      'validphys.results',
                      'validphys.covmats',
                      'n3fit.n3fit_checks_provider'
]

SETUPFIT_DEFAULTS = dict(
    use_cuts = 'internal',
)



log = logging.getLogger(__name__)

RUNCARD_COPY_FILENAME = "filter.yml"
FILTER_OUTPUT_FOLDER = "filter"
TABLE_OUTPUT_FOLDER = "tables"
MD5_FILENAME = "md5"
INPUT_FOLDER = "input"


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
            if not self.config_yml.is_file():
                raise SetupFitError("Invalid runcard. Must be a file.")

        # check if results folder exists
        self.output_path = pathlib.Path(self.output_path).absolute()

        if self.output_path.is_dir():
            log.warning(f"Output folder exists: {self.output_path} Overwriting contents")
        else:
            if not re.fullmatch(r'[\w\-]+', self.output_path.name):
                raise SetupFitError("Invalid output folder name. Must be alphanumeric.")
            try:
                self.output_path.mkdir()
            except OSError as e:
                raise EnvironmentError_(e) from e

        try:
            shutil.copy2(self.config_yml, self.output_path / RUNCARD_COPY_FILENAME)
        except shutil.SameFileError:
            pass
        except Exception as e:
            raise EnvironmentError_(e) from e

        # create output folder
        self.filter_path = self.output_path / FILTER_OUTPUT_FOLDER
        self.filter_path.mkdir(exist_ok=True)
        self.table_folder = self.output_path / TABLE_OUTPUT_FOLDER
        self.table_folder.mkdir(exist_ok=True)
        # put lockfile input inside of filter output
        self.input_folder = self.filter_path / INPUT_FOLDER
        self.input_folder.mkdir(exist_ok=True)

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

        if file_content.get('closuretest') is not None:
            filter_action = 'datacuts::closuretest::theory::fitting filter'
            check_n3fit_action = 'datacuts::theory::closuretest::fitting n3fit_checks_action'
        else:
            filter_action = 'datacuts::theory::fitting filter'
            check_n3fit_action = 'datacuts::theory::fitting n3fit_checks_action'
        
        SETUPFIT_FIXED_CONFIG['actions_'] += [check_n3fit_action, filter_action]
        
        if file_content.get('theorycovmatconfig') is not None:
            SETUPFIT_FIXED_CONFIG['actions_'].append(
                'datacuts::theory::theorycovmatconfig nnfit_theory_covmat')
        for k,v in SETUPFIT_DEFAULTS.items():
            file_content.setdefault(k, v)
        file_content.update(SETUPFIT_FIXED_CONFIG)
        return cls(file_content, *args, **kwargs)


class SetupFitApp(App):
    """The class which parsers and perform the filtering"""
    environment_class = SetupFitEnvironment
    config_class = SetupFitConfig

    def __init__(self):
        super(SetupFitApp, self).__init__(name='setup-fit',
                                          providers=SETUPFIT_PROVIDERS)

    @property
    def argparser(self):
        parser = super().argparser
        parser.add_argument('-o','--output',
                        help="Output folder and name of the fit",
                        default=None)
        parser.add_argument("--legacy",
                            help="Filter an old nnfit runcard by skipping n3fit specific checks",
                            action='store_true')
        return parser

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        if args['output'] is None:
            args['output'] = pathlib.Path(args['config_yml']).stem
        return args

    def run(self):
        try:
            # set folder output name
            self.environment.config_yml = pathlib.Path(self.args['config_yml']).absolute()
            self.environment.legacy = self.args["legacy"]

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
