"""
This module contains the CLI for evolven3fit
"""
import logging
import pathlib
import numpy as np
from argparse import ArgumentParser

from evolven3fit import evolven3fit_API
from validphys.core import CommonDataMetadata

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser(
        description="evolven3fit - a script with tools to evolve PDF fits"
    )
    subparsers = parser.add_subparsers(title="actions", dest="actions")
    evolvefit_parser = subparsers.add_parser(
        "evolve",
        help="Evolves the fitted PDFs. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. If a path is given for the dump option, the eko will be dumped in that path after the computation. If a path is given for the load option, the eko to be used for the evolution will be loaded from that path. The two options are mutually exclusive.",
    )
    eko_parser = subparsers.add_parser(
        "produce_eko",
        help="Produce the eko for the specified theory. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. The x_grid starts at x_grid_ini and ends at 1.0 and contains the provided number of points. The eko will be dumped in the provided path.",
    )
    evolvefit_parser.add_argument(
        "configuration_folder",
        help="Path to the folder containing the (pre-DGLAP) fit result",
    )
    parser.add_argument(
        "-q",
        "--q_fin",
        type=float,
        default=None,
        help="Final q-value of the evolution",
    )
    parser.add_argument(
        "-p",
        "--q_points",
        type=int,
        default=None,
        help="Number of q points for the evolution",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        default=1,
        help="Number of cores to be used",
    )
    evolvefit_parser.add_argument(
        "-l",
        "--load",
        type=pathlib.Path,
        default=None,
        help="Path of the EKO to be loaded",
    )
    evolvefit_parser.add_argument(
        "-d",
        "--dump",
        type=pathlib.Path,
        default=None,
        help="Path where the EKO is dumped (optional)",
    )
    eko_parser.add_argument(
        "theoryID", type=int, help="ID of the theory used to produce the eko"
    )
    eko_parser.add_argument(
        "dump",
        type=pathlib.Path,
        help="Path where the EKO is dumped",
    )
    eko_parser.add_argument(
        "x_grid_ini",
        default=1.0e-7,
        type=float,
        help="Starting point of the x-grid",
    )
    eko_parser.add_argument(
        "x_grid_points",
        default=196,
        type=int,
        help="Number of points of the x-grid",
    )
    args = parser.parse_args()
    op_card_info = {
        "ev_op_max_order": 10,
        "ev_op_iterations": 1,
        "n_integration_cores": args.n_cores,
        "backward_inversion": "expanded",
    }
    t_card_info = {}
    if args.actions == "evolve":
        cli_evolven3fit(
            args.configuration_folder,
            args.q_fin,
            args.q_points,
            op_card_info,
            t_card_info,
            args.dump,
            args.load,
        )
    elif args.actions == "produce_eko":
        x_grid = np.geomspace(args.x_grid_ini, 1.0, args.x_grid_points)
        tcard, opcard = evolven3fit_API.construct_eko_cards(
            args.theoryID, op_card_info, t_card_info, args.q_fin, args.q_points, x_grid
        )
        eko_op = evolven3fit_API.construct_eko_for_fit(tcard, opcard, args.dump)


def cli_evolven3fit(
    configuration_folder, q_fin, q_points, op_card_info, t_card_info, dump, load
):
    """Evolves the fitted PDFs.

    The q_grid starts at the Q0 given by the theory but
    the last point is q_fin and its number of
    points can be specified by q_points. If just one of the
    two is not specified by the user, the default grid
    will be used.

    If a path is given for the dump option, the eko
    will be dumped in that path after the computation.

    If a path is given for the load option, the eko
    to be used for the evolution will be loaded from that
    path.

    The two options are mutually exclusive.
    """
    return evolven3fit_API.evolve_fit(
        configuration_folder, q_fin, q_points, op_card_info, t_card_info, load, dump
    )


if __name__ == "__main__":
    main()
