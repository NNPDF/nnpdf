# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:57:38 2016

@author: Zahari Kassabov
"""
#TODO: Move most of this to reportengine
import argparse
import logging
import pathlib
import sys
import contextlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


from reportengine.resourcebuilder import ResourceBuilder, ResourceError

from validphys import providers
from validphys.config import Config, ConfigError, Environment

log = logging.getLogger(__name__)

def format_rich_error(e):
    with contextlib.redirect_stdout(sys.stderr):
        log.error("Bad configuration encountered:")
        print(e)


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

    parser.add_argument('--style',
                        help='matplotlib style file to override the built-in one.',
                        default=None)


    parser.add_argument('--formats', nargs='+', help="formats of the output figures",
                        default=('pdf',))

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
                              this_folder = pathlib.Path(__file__).parent,
                              figure_formats=args.formats
                              )
    environment.init_output()

    try:
        with open(args.config_yml) as f:
            try:
                c = Config.from_yaml(f, environment=environment)
            except ConfigError as e:
                format_rich_error(e)
                sys.exit(1)

    except OSError as e:
        log.error("Could not open configuration file: %s" % e)
        sys.exit(1)

    if args.style:
        try:
            plt.style.use(args.style)
        except Exception as e:
            log.error("There was a problem reading the supplied style: %s" %e,
                  file=sys.stderr)
            sys.exit(1)
    else:
        plt.style.use(str(environment.this_folder / 'small.mplstyle'))


    try:
        actions = c.parse_actions_(c['actions_'])
    except ConfigError as e:
        format_rich_error(e)
        sys.exit(1)
    except KeyError as e:
        log.error("A key 'actions_' is needed in the top level of the config file.")
        sys.exit(1)

    rb = ResourceBuilder(c, providers, actions, environment=environment)

    #environment.output_path.mkdir(exist_ok = True)

    try:
        rb.resolve_targets()
    except ConfigError as e:
        format_rich_error(e)
        sys.exit(1)
    except ResourceError as e:
        with contextlib.redirect_stdout(sys.stderr):
            log.error("Cannot process a resource:")
            print(e)
        sys.exit(1)
    rb.execute_sequential()
