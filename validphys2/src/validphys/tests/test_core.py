"""
Test core functionality
"""
import pytest
from validphys import core
from validphys.tests.conftest import PDF, HESSIAN_PDF

@pytest.mark.parametrize("pdf_name", [PDF, HESSIAN_PDF])
def test_pdf(pdf_name):
    """Check that the given PDF and their relevant attributes can be read
    And that they don't have crazy values
    """
    pdf = core.PDF(pdf_name)
    _ = pdf.q_min
    _ = pdf.stats_class
    assert pdf.isinstalled
    error_type = pdf.error_type
    if error_type == "replicas":
        assert pdf.get_members() == (len(pdf)-1)
        assert pdf.error_conf_level is None
    else:
        assert pdf.get_members() == len(pdf)
        assert isinstance(pdf.error_conf_level, (int, float))
    assert pdf.name == pdf._plotname == pdf_name == str(pdf)
