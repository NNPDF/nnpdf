"""
This module contains the CLI for evolven3fit
"""

from argparse import ArgumentParser
import logging
import pathlib
import sys

from evolven3fit import cli, eko_utils, evolve, utils
import numpy as np

from eko.runner.managed import solve
from n3fit.io.writer import XGRID
from nnpdf_data import theory_cards
from nnpdf_data.theorydbutils import fetch_theory
from validphys.loader import FallbackLoader, Loader

_logger = logging.getLogger(__name__)


def _add_production_options(parser):
    """Adds the production options that are shared between the normal and the QED EKOs to their parsers"""
    parser.add_argument("theoryID", type=int, help="ID of the theory used to produce the eko")
    parser.add_argument("dump", type=pathlib.Path, help="Path where the EKO is dumped")
    parser.add_argument(
        "-i", "--x-grid-ini", default=None, type=float, help="Starting point of the x-grid"
    )
    parser.add_argument(
        "-p", "--x-grid-points", default=None, type=int, help="Number of points of the x-grid"
    )
    parser.add_argument("-n", "--n-cores", type=int, default=1, help="Number of cores to be used")
    parser.add_argument('--use_polarized', action='store_true', help="Use polarized evolution")
    parser.add_argument(
        "-e",
        "--ev-op-iterations",
        type=int,
        default=None,
        help="ev_op_iterations for the EXA theory. Overrides the settings given in the theory card.",
    )


def construct_eko_parser(subparsers):
    parser = subparsers.add_parser(
        "produce_eko",
        help=(
            """Produce the eko for the specified theory.
            The q_grid starts at the Q0 given by the theory but the last point is q_fin
            and its number of points can be specified by q_points.
            The x_grid starts at x_grid_ini and ends at 1.0 and contains the
            provided number of points. The eko will be dumped in the provided path."""
        ),
    )
    _add_production_options(parser)
    parser.add_argument(
        "--legacy40",
        action="store_true",
        help="Use evolution grid used in NNPDF4.0 (for reproducibility)",
    )
    return parser


def construct_eko_photon_parser(subparsers):
    parser = subparsers.add_parser(
        "produce_eko_photon",
        help=(
            """Produce the eko_photon for the specified theory.
            The q_gamma will be the provided one.
            The x_grid starts at x_grid_ini and ends at 1.0 and contains the
            provided number of points. The eko will be dumped in the provided path."""
        ),
    )
    _add_production_options(parser)
    parser.add_argument(
        "-g", "--q-gamma", default=100, type=float, help="Scale at which the photon is generated"
    )
    return parser


def construct_evolven3fit_parser(subparsers):
    parser = subparsers.add_parser(
        "evolve",
        help="Evolves the fitted PDFs. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. If a path is given for the dump option, the eko will be dumped in that path after the computation. If a path is given for the load option, the eko to be used for the evolution will be loaded from that path. The two options are mutually exclusive.",
    )
    parser.add_argument(
        "fit_folder", help="Path to the folder containing the (pre-DGLAP) fit result"
    )
    parser.add_argument(
        "-l", "--load", type=pathlib.Path, default=None, help="Path of the EKO to be loaded"
    )
    parser.add_argument(
        "-d", "--dump", type=pathlib.Path, default=None, help="Path where the EKO is dumped"
    )
    parser.add_argument(
        "--hessian_fit", type=bool, default=False, help="Specify if the fit is hessian (default is False)"
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force the evolution to be done even if it has already been done",
    )
    return parser


def evolven3fit_new():
    _logger.warning("`evolven3fit_new` is deprecated. Please use `evolven3fit` instead.")
    main()


def main():
    parser = ArgumentParser(description="evolven3fit - a script with tools to evolve PDF fits")
    parser.add_argument(
        "-q", "--q-fin", type=float, default=None, help="Final q-value of the evolution"
    )
    parser.add_argument(
        "-p", "--q-points", type=int, default=None, help="Number of q points for the evolution"
    )
    parser.add_argument("--no-net", action="store_true", help="Emulates validphys' --no-net")

    subparsers = parser.add_subparsers(title="actions", dest="actions")
    construct_eko_parser(subparsers)
    construct_eko_photon_parser(subparsers)
    construct_evolven3fit_parser(subparsers)

    args = parser.parse_args()

    op_card_info = {}
    # Here we do not allow any modification of the theory card, for the moment.
    theory_card_info = {}

    if args.no_net:
        loader = Loader()
    else:
        loader = FallbackLoader()

    if args.actions == "evolve":

        fit_folder = pathlib.Path(args.fit_folder)
        if args.load is None:
            # no path provided to load the eko
            utils.check_filter(fit_folder)
            theoryID = utils.get_theoryID_from_runcard(fit_folder)
            _logger.info(f"Loading eko from theory {theoryID}")
            try:
                # try to load eko from the theory
                eko_path = loader.check_eko(theoryID)
            except:
                # if the theory is not there and a path is given 
                # for the dump option, the eko will be dumped in that path after the computation.
                eko_path = evolve.create_eko (
                    fit_folder, 
                    args.args.q_fin, 
                    args.q_points, 
                    op_card_info, 
                    theory_card_info, 
                    args.dump,
                )

        else:
            # If a path is given for the load option, the eko
            # to be used for the evolution will be loaded from that path.
            eko_path = args.load

        cli.cli_evolven3fit(
            fit_folder,
            args.force,
            eko_path,
            args.hessian_fit
        )
    else:
        # If we are in the business of producing an eko, do some checks before starting:
        # 1. load the nnpdf theory early to check for inconsistent options and theory problems
        nnpdf_theory = fetch_theory(theory_cards, args.theoryID)
        if nnpdf_theory.get("ModEv") == "TRN" and args.ev_op_iterations is not None:
            raise ValueError("ev_op_iterations is not accepted with ModEv=TRN solution")

        stdout_log = logging.StreamHandler(sys.stdout)
        stdout_log.setLevel(evolve.LOGGING_SETTINGS["level"])
        stdout_log.setFormatter(evolve.LOGGING_SETTINGS["formatter"])
        for logger_ in (_logger, *[logging.getLogger("eko")]):
            logger_.handlers = []
            logger_.setLevel(evolve.LOGGING_SETTINGS["level"])
            logger_.addHandler(stdout_log)
        if args.x_grid_ini is None:
            if args.x_grid_points is None:
                x_grid = XGRID
            else:
                raise ValueError(
                    "x_grid_ini and x_grid_points must be specified either both or none of them"
                )
        elif args.x_grid_points is None:
            raise ValueError(
                "x_grid_ini and x_grid_points must be specified either both or none of them"
            )
        else:
            x_grid = np.geomspace(args.x_grid_ini, 1.0, args.x_grid_points)

        # Prepare the op card config
        op_card_info = {
            "configs": {"n_integration_cores": args.n_cores, "polarized": args.use_polarized}
        }
        if args.ev_op_iterations is not None:
            op_card_info["configs"]["ev_op_iterations"] = args.ev_op_iterations

        if args.actions == "produce_eko":
            tcard, opcard = eko_utils.construct_eko_cards(
                nnpdf_theory,
                args.q_fin,
                args.q_points,
                x_grid,
                op_card_info,
                theory_card_info,
                args.legacy40,
            )
        elif args.actions == "produce_eko_photon":
            tcard, opcard = eko_utils.construct_eko_photon_cards(
                nnpdf_theory, args.q_fin, x_grid, args.q_gamma, op_card_info, theory_card_info
            )
        solve(tcard, opcard, args.dump)


if __name__ == "__main__":
    main()
