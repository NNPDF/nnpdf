# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:57:38 2016

@author: Zahari Kassabov
"""
import argparse
import logging
import pathlib
import sys

import matplotlib.pyplot as plt

from reportengine.resourcebuilder import ResourceBuilder, ResourceError

from validphys import providers
from validphys.config import Config, ConfigError, Environment


def main():
    #TODO: Oberhaul this to use reportengine properly

    parser = argparse.ArgumentParser(
             description = "Validphys developer preview",
             )

    parser.add_argument('config_yml',
                        help = "path to the configuration file")

    parser.add_argument('-o','--output', help="output folder where to "
                                         "store resulting plots and tables",
                        default='output')

    loglevel = parser.add_mutually_exclusive_group()

    loglevel.add_argument('-q','--quiet', help="supress INFO messages and C output",
                        action='store_true')

    loglevel.add_argument('-d', '--debug', help = "show debug info",
                          action='store_true')

    parser.add_argument('-p','--datapath', help="path where the NNPDF "
                        "data is located",
                        default='../nnpdfcpp/data')

    parser.add_argument('--resultspath', help="path where the fit results "
                          "are located. Calculated from 'datapath' by default",
                         )

    args = parser.parse_args()

    if args.quiet:
        level = logging.WARN
    elif args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(format='%(levelname)s: %(message)s', level=level)

    if not args.resultspath:
        args.resultspath = pathlib.Path(args.datapath).parent / 'nnpdfbuild' / 'results'

    environment = Environment(data_path = args.datapath,
                              results_path = args.resultspath,
                              output_path = args.output,
                              this_folder = pathlib.Path(__file__).parent
                              )
    environment.init_output()

    try:
        with open(args.config_yml) as f:
            try:
                c = Config.from_yaml(f, environment=environment)
            except ConfigError as e:
                print("Bad configuration encountered:")
                print(e)
                print(e.alternatives_text())
                sys.exit(1)

    except OSError as e:
        print("Could not open configuration file: %s" % e)
        sys.exit(1)


    actions = c.parse_actions_(c['actions_'])

    rb = ResourceBuilder(c, providers, actions)
    rb.rootns['environment'] = environment

    #environment.output_path.mkdir(exist_ok = True)
    plt.style.use(str(environment.this_folder / 'small.mplstyle'))
    try:
        rb.resolve_targets()
    except ConfigError as e:
        print("Bad configuration encountered:")
        print(str(e.args[0]))
        print(e.alternatives_text())
        sys.exit(1)
    except ResourceError as e:
        print("Cannot process a resource:")
        print(e)
    rb.execute_sequential()
