"""
This module contains the CLI for evolven3fit_new
"""
import logging
import pathlib
from argparse import ArgumentParser
import sys
import numpy as np

from evolven3fit_new import eko_utils, evolve, cli

_logger = logging.getLogger(__name__)


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
    parser.add_argument(
        "theoryID", type=int, help="ID of the theory used to produce the eko"
    )
    parser.add_argument(
        "dump",
        type=pathlib.Path,
        help="Path where the EKO is dumped",
    )
    parser.add_argument(
        "-i",
        "--x-grid-ini",
        default=None,
        type=float,
        help="Starting point of the x-grid",
    )
    parser.add_argument(
        "-p",
        "--x-grid-points",
        default=None,
        type=int,
        help="Number of points of the x-grid",
    )
    return parser


def construct_evolven3fit_parser(subparsers):
    parser = subparsers.add_parser(
        "evolve",
        help="Evolves the fitted PDFs. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. If a path is given for the dump option, the eko will be dumped in that path after the computation. If a path is given for the load option, the eko to be used for the evolution will be loaded from that path. The two options are mutually exclusive.",
    )
    parser.add_argument(
        "configuration_folder",
        help="Path to the folder containing the (pre-DGLAP) fit result",
    )
    parser.add_argument(
        "-l",
        "--load",
        type=pathlib.Path,
        default=None,
        help="Path of the EKO to be loaded",
    )
    parser.add_argument(
        "-d",
        "--dump",
        type=pathlib.Path,
        default=None,
        help="Path where the EKO is dumped",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force the evolution to be done even if it has already been done",
    )
    return parser


def main():
    parser = ArgumentParser(
        description="evolven3fit_new - a script with tools to evolve PDF fits"
    )
    parser.add_argument(
        "-q",
        "--q-fin",
        type=float,
        default=None,
        help="Final q-value of the evolution",
    )
    parser.add_argument(
        "-p",
        "--q-points",
        type=int,
        default=None,
        help="Number of q points for the evolution",
    )
    parser.add_argument(
        "-n",
        "--n-cores",
        type=int,
        default=1,
        help="Number of cores to be used",
    )
    subparsers = parser.add_subparsers(title="actions", dest="actions")
    eko_parser = construct_eko_parser(subparsers)
    evolven3fit_parser = construct_evolven3fit_parser(subparsers)

    args = parser.parse_args()
    op_card_info = {
        "ev_op_max_order": 10,
        "ev_op_iterations": 1,
        "n_integration_cores": args.n_cores,
        "backward_inversion": "expanded",
    }
    theory_card_info = {}
    if args.actions == "evolve":
        cli.cli_evolven3fit_new(
            args.configuration_folder,
            args.q_fin,
            args.q_points,
            op_card_info,
            theory_card_info,
            args.dump,
            args.load,
            args.force,
        )
    elif args.actions == "produce_eko":
        stdout_log = logging.StreamHandler(sys.stdout)
        stdout_log.setLevel(evolve.LOGGING_SETTINGS["level"])
        stdout_log.setFormatter(evolve.LOGGING_SETTINGS["formatter"])
        for logger_ in (_logger, *[logging.getLogger("eko")]):
            logger_.handlers = []
            logger_.setLevel(evolve.LOGGING_SETTINGS["level"])
            logger_.addHandler(stdout_log)
        if args.x_grid_ini is None:
            if args.x_grid_points is None:
                x_grid = cli.XGRID
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
        tcard, opcard = eko_utils.construct_eko_cards(
            args.theoryID, args.q_fin, args.q_points, x_grid, op_card_info, theory_card_info
        )
        eko_op = eko_utils.construct_eko_for_fit(tcard, opcard, args.dump)


if __name__ == "__main__":
    main()
