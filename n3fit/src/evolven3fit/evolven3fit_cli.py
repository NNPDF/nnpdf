"""
This module contains the CLI for evolven3fit
"""
import logging
import pathlib
from argparse import ArgumentParser

from evolven3fit import evolven3fit_API

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser(description="evolven3fit - a script to evolve PDF fits")
    parser.add_argument(
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
    args = parser.parse_args()
    cli_evolven3fit(
        args.configuration_folder,
        args.q_fin,
        args.q_points,
        10,
        1,
        args.n_cores,
        "expanded",
        args.dump,
        args.load,
    )


def cli_evolven3fit(
    configuration_folder, q_fin, q_points, ev_ord, ev_it, n_cores, back, dump, load
):
    """Evolves the fitted PDFs.

    The q_grid starts at the Q0 given by the theory but
    the last point is q_fin and its number of
    points can be specified by q_points. If just one of the
    two is not specified by the user, the default grid
    will be used.

    If a path is given for the dump option, the eko
    will be dumped in that path after the computation.

    If a path is given for the dump option, the eko
    to be used for the evolution will be loaded from that
    path.

    The two options are mutually exclusive.
    """
    op_card_info = {
        "ev_op_max_order": ev_ord,
        "ev_op_iterations": ev_it,
        "n_integration_cores": n_cores,
        "backward_inversion": back,
    }
    t_card_info = {}
    import ipdb

    ipdb.set_trace()
    return evolven3fit_API.evolve_fit(
        configuration_folder, q_fin, q_points, op_card_info, t_card_info, load, dump
    )
