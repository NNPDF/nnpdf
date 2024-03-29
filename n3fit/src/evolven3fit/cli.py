from . import evolve, utils


def cli_evolven3fit(
    configuration_folder, q_fin, q_points, op_card_info, theory_card_info, dump, load, force
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
    utils.check_is_a_fit(configuration_folder)
    return evolve.evolve_fit(
        configuration_folder, q_fin, q_points, op_card_info, theory_card_info, force, load, dump
    )
