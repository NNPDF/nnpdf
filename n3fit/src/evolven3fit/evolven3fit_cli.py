"""
This module contains the CLI for evolven3fit
"""

import click

from . import evolven3fit_API


@click.group()
def cli():
    pass


@cli.command("evolve")
@click.argument("configuration_folder", nargs=1)
def cli_evolven3fit(configuration_folder):
    """Evolves the fitted PDFs"""
    return evolven3fit_API.evolve_fit(configuration_folder)
