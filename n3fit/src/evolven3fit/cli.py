from . import evolve, utils


def cli_evolven3fit(configuration_folder, force, load, hessian):
    """Evolves the fitted PDFs."""

    utils.check_nnfit_folder(configuration_folder)
    return evolve.evolve_fit(configuration_folder, force, load, hessian)
