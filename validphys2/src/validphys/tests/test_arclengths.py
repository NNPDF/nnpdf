from validphys.api import API
from validphys.tableloader import sane_load
from validphys.tests.test_regressions import make_table_comp


@make_table_comp(sane_load)
def test_arclength_mc(mc_pdf_config):
    """Integration test that ``arc_length_table`` action matches expected value
    when using a MC ``PDF``.

    This checks that the underlying ``arc_lengths`` action successfully
    runs and produces expected numbers.

    """
    table = API.arc_length_table(**mc_pdf_config)
    return table


@make_table_comp(sane_load)
def test_arclength_hessian(hessian_pdf_config):
    """Like ``test_arclength_mc`` except uses a hessian ``PDF``. Checks that
    the underlying action produces correct result when errortype is hessian.

    """
    table = API.arc_length_table(**hessian_pdf_config)
    return table
