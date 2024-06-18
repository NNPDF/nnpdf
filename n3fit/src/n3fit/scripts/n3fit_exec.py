#!/usr/bin/env python
"""
    n3fit - performs fit using ml external frameworks
"""

import argparse
import logging
import pathlib
import re
import shutil
import sys
import warnings

from reportengine import colors
from reportengine.compat import yaml
from reportengine.namespaces import NSList
from validphys.app import App
from validphys.config import Config, ConfigError, Environment, EnvironmentError_
from validphys.core import FitSpec

N3FIT_FIXED_CONFIG = dict(use_cuts='internal', use_t0=True, actions_=[])

FIT_NAMESPACE = "datacuts::theory::fitting "
CLOSURE_NAMESPACE = "datacuts::theory::closuretest::fitting "

N3FIT_PROVIDERS = [
    "n3fit.performfit",
    "n3fit.n3fit_checks_provider",
    "validphys.results",
    "validphys.n3fit_data",
    "validphys.pseudodata",
    "validphys.covmats",
    "validphys.commondata",
]

log = logging.getLogger(__name__)

RUNCARD_COPY_FILENAME = "filter.yml"
INPUT_FOLDER = "input"
TAB_FOLDER = "tables"


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
        if not (self.output_path / "nnfit").is_dir():
            if not re.fullmatch(r"[\w.\-]+", self.output_path.name):
                raise N3FitError("Invalid output folder name. Must be alphanumeric.")
            try:
                self.output_path.mkdir(exist_ok=True)
                (self.output_path / "nnfit").mkdir(exist_ok=True)
            except OSError as e:
                raise EnvironmentError_(e) from e

            try:
                shutil.copy2(self.config_yml, self.output_path / RUNCARD_COPY_FILENAME)
            except shutil.SameFileError:
                pass

        # create output folder for the fit
        self.replica_path = self.output_path / "nnfit"
        for replica in self.replicas:
            path = self.replica_path / "replica_{0}".format(replica)
            log.info("Creating replica output folder in {0}".format(path))
            try:
                path.mkdir(exist_ok=True)
            except OSError as e:
                raise EnvironmentError_(e) from e

        # place tables in last replica path, folder already exists.
        self.table_folder = path

        # make lockfile input inside of replica folder
        # avoid conflict with setupfit
        self.input_folder = self.replica_path / INPUT_FOLDER
        self.input_folder.mkdir(exist_ok=True)

    @classmethod
    def ns_dump_description(cls):
        return {
            "replicas": "The MC replica number/s",
            "replica_path": "The replica output path",
            "output_path": "The runcard name",
            "hyperopt": "The hyperopt flag",
            **super().ns_dump_description(),
        }


class N3FitConfig(Config):
    """Specialization for yaml parsing"""

    @classmethod
    def from_yaml(cls, o, *args, **kwargs):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", yaml.error.MantissaNoDotYAML1_1Warning)
                # We need to specify the older version 1.1 to support the
                # older configuration files, which liked to use on/off for
                # booleans.
                # The floating point parsing yields warnings everywhere, which
                # we suppress.
                file_content = yaml.safe_load(o, version="1.1")
        except yaml.error.YAMLError as e:
            raise ConfigError(f"Failed to parse yaml file: {e}")
        if not isinstance(file_content, dict):
            raise ConfigError(
                f"Expecting input runcard to be a mapping, not '{type(file_content)}'."
            )

        if file_content.get('closuretest') is not None:
            namespace = CLOSURE_NAMESPACE
        else:
            namespace = FIT_NAMESPACE

        N3FIT_FIXED_CONFIG['actions_'].append(namespace + "performfit")

        if fps := file_content["fitting"].get("savepseudodata", True):
            if fps != True:
                raise TypeError(f"fitting::savepseudodata is neither True nor False ({fps})")
            if len(kwargs["environment"].replicas) != 1:
                raise ConfigError(
                    "Cannot request that multiple replicas are fitted and that "
                    "pseudodata is saved. Either set `fitting::savepseudodata` "
                    "to `false` or fit replicas one at a time."
                )
            # take same namespace configuration on the pseudodata_table action.
            training_action = namespace + "training_pseudodata"
            validation_action = namespace + "validation_pseudodata"

            N3FIT_FIXED_CONFIG['actions_'].extend((training_action, validation_action))

        if thconfig := file_content.get('fiatlux'):
            N3FIT_FIXED_CONFIG['fiatlux'] = thconfig
        else:
            N3FIT_FIXED_CONFIG['fiatlux'] = None

        if thconfig := file_content.get('positivity_bound'):
            N3FIT_FIXED_CONFIG['positivity_bound'] = thconfig
        else:
            N3FIT_FIXED_CONFIG['positivity_bound'] = None

        # Theorycovmat flags and defaults
        N3FIT_FIXED_CONFIG['theory_covmat_flag'] = False
        N3FIT_FIXED_CONFIG['use_thcovmat_in_fitting'] = False
        N3FIT_FIXED_CONFIG['use_thcovmat_in_sampling'] = False
        if (thconfig := file_content.get('theorycovmatconfig')) is not None:
            N3FIT_FIXED_CONFIG['use_thcovmat_in_fitting'] = thconfig.get(
                'use_thcovmat_in_fitting', True
            )
            N3FIT_FIXED_CONFIG['use_thcovmat_in_sampling'] = thconfig.get(
                'use_thcovmat_in_sampling', True
            )
            if (
                N3FIT_FIXED_CONFIG['use_thcovmat_in_sampling']
                or N3FIT_FIXED_CONFIG['use_thcovmat_in_fitting']
            ):
                N3FIT_FIXED_CONFIG['theory_covmat_flag'] = True
            N3FIT_FIXED_CONFIG['use_user_uncertainties'] = thconfig.get(
                'use_user_uncertainties', False
            )
            N3FIT_FIXED_CONFIG['use_scalevar_uncertainties'] = thconfig.get(
                'use_scalevar_uncertainties', True
            )
        # Sampling flags
        if (sam_t0 := file_content.get('sampling')) is not None:
            N3FIT_FIXED_CONFIG['separate_multiplicative'] = sam_t0.get(
                'separate_multiplicative', False
            )
        # Fitting flag
        file_content.update(N3FIT_FIXED_CONFIG)
        return cls(file_content, *args, **kwargs)

    def produce_fit(self):
        """Produces a FitSpec which points at the current fit, to load
        fit_commondata from.
        """
        fitpath = self.environment.output_path
        return FitSpec(fitpath.name, fitpath)

    def parse_fakedata(self, fakedata: bool):
        """Parses the `fakedata` key from the closuretest namespace, if True then
        use generated closure test in fit
        """
        if fakedata:
            log.warning("using filtered closure data")
            if not (self.environment.output_path / 'filter').is_dir():
                raise ConfigError(
                    "Could not find filter result at "
                    f"{self.environment.output_path/'filter'} "
                    "to load commondata from. Did you run filter on the "
                    "runcard?"
                )
        return fakedata

    def produce_use_fitcommondata(self, fakedata):
        """Produces the `use_fitcommondata` key from the `fakedata` key in
        `closuretest` namespace
        """
        return fakedata

    def produce_kfold_parameters(self, kfold=None, hyperopt=None):
        """Return None even if there are kfolds in the runcard if the hyperopt flag is not active"""
        if hyperopt is not None:
            return kfold
        return None

    def produce_kpartitions(self, kfold_parameters):
        if kfold_parameters:
            partitions = kfold_parameters["partitions"]
            # Note that one of the partitions could be empty ([]) or, by yaml usual notation, None
            for partition in partitions:
                if partition["datasets"] is None:
                    partition["datasets"] = []
            return partitions
        return None

    def produce_hyperscanner(self, parameters, hyperscan_config=None, hyperopt=None):
        """For a hyperparameter scan to be run, a hyperscanner must be
        constructed from the original hyperscan_config"""
        # This needs to be imported here because it needs Tensorflow and n3fit
        from n3fit.hyper_optimization.hyper_scan import HyperScanner

        if hyperscan_config is None or hyperopt is None:
            return None
        if hyperopt and self.environment.restart:
            hyperscan_config.update({'restart': 'true'})
        if hyperopt and self.environment.parallel_hyperopt:
            hyperscan_config.update({'parallel': 'true'})
            hyperscan_config.update(
                {
                    'db_host': self.environment.db_host,
                    'db_port': self.environment.db_port,
                    'db_name': self.environment.db_name,
                    'output_path': self.environment.output_path.name,
                    'num_mongo_workers': self.environment.num_mongo_workers,
                }
            )
        return HyperScanner(parameters, hyperscan_config)


class N3FitApp(App):
    """The class which parsers and performs the fit"""

    environment_class = N3FitEnvironment
    config_class = N3FitConfig

    def __init__(self):
        super(N3FitApp, self).__init__(name="n3fit", providers=N3FIT_PROVIDERS)

    @property
    def argparser(self):
        parser = super().argparser
        parser.add_argument(
            "-o", "--output", help="Output folder and name of the fit", default=None
        )

        def check_positive(value):
            ivalue = int(value)
            if ivalue <= 0:
                raise argparse.ArgumentTypeError("%s is an invalid positive int value." % value)
            return ivalue

        parser.add_argument("--hyperopt", help="Enable hyperopt scan", default=None, type=int)
        parser.add_argument("--restart", help="Enable hyperopt restarts", action="store_true")
        parser.add_argument(
            "--parallel-hyperopt",
            help="Enable hyperopt run in parallel with MongoDB",
            action="store_true",
        )
        parser.add_argument("--db-host", help="MongoDB host", default="localhost")
        parser.add_argument("--db-port", help="MongoDB port", default=27017)
        parser.add_argument("--db-name", help="MongoDB dataset name", default="hyperopt-db")
        parser.add_argument(
            "--num-mongo-workers",
            help="Number of mongo workers to be launched simultaneously",
            type=check_positive,
            default=1,
        )
        parser.add_argument("replica", help="MC replica number", type=check_positive)
        parser.add_argument(
            "-r",
            "--replica_range",
            help="End of the range of replicas to compute",
            type=check_positive,
        )
        return parser

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)

        # Validate dependencies related to the --hyperopt argument
        if args["hyperopt"] is None:
            if args["restart"]:
                raise argparse.ArgumentError(
                    None, "The --restart option requires --hyperopt to be set."
                )
            if args["parallel_hyperopt"]:
                raise argparse.ArgumentError(
                    None, "The --parallel-hyperopt option requires --hyperopt to be set."
                )

        if args["output"] is None:
            args["output"] = pathlib.Path(args["config_yml"]).stem
        return args

    def run(self):
        try:
            self.environment.config_yml = pathlib.Path(self.args["config_yml"]).absolute()
            replica = self.args["replica"]
            if self.args["replica_range"]:
                replicas = list(range(replica, self.args["replica_range"] + 1))
            else:
                replicas = [replica]
            self.environment.replicas = NSList(replicas, nskey="replica")
            self.environment.hyperopt = self.args["hyperopt"]
            self.environment.restart = self.args["restart"]
            self.environment.parallel_hyperopt = self.args["parallel_hyperopt"]
            self.environment.db_host = self.args["db_host"]
            self.environment.db_port = self.args["db_port"]
            self.environment.db_name = self.args["db_name"]
            self.environment.num_mongo_workers = self.args["num_mongo_workers"]
            super().run()
        except N3FitError as e:
            log.error(f"Error in n3fit:\n{e}")
            sys.exit(1)
        except Exception as e:
            log.critical(f"Bug in n3fit ocurred. Please report it.")
            print(colors.color_exception(e.__class__, e, e.__traceback__), file=sys.stderr)
            sys.exit(1)


def main():
    a = N3FitApp()
    a.main()


if __name__ == "__main__":
    main()
