"""
This module contains the CLI for evolven3fit
"""

import pathlib
import click

from . import evolven3fit_API


@click.group()
def cli():
    pass


@cli.command("evolve")
@click.argument("configuration_folder", nargs=1)
@click.option("-o", "--ev_ord", type=int, default = 10)
@click.option("-i", "--ev_it", type=int, default = 1)
@click.option("-n", "--n_cores", type=int, default = 1)
@click.option("-b", "--back", type=str, default = "expanded")
@click.option("-d", "--dump", type=pathlib.Path, default = None)
@click.option("-l", "--load", type=pathlib.Path, default = None)
def cli_evolven3fit(configuration_folder, ev_ord, ev_it, n_cores, back, dump, load):
    """Evolves the fitted PDFs.

        If a path is given for the dump option, the eko 
        will be dumped in that path after the computation.

        If a path is given for the dump option, the eko
        to be used for the evolution will be loaded from that
        path.

        The two options are mutually exclusive.
    """
    op_card_info = {"ev_op_max_order": ev_ord, "ev_op_iterations": ev_it, "n_integration_cores": n_cores, "backward_inversion": back}
    t_card_info = {}
    return evolven3fit_API.evolve_fit(configuration_folder, op_card_info, t_card_info, load, dump)
