import lhapdf

from validphys.api import API
from validphys.tests.conftest import PDF, temp_lhapdf_path


def test_alpha_s_bundle_pdf(tmp_path):
    API.alpha_s_bundle_pdf(pdf=PDF, pdfs=[PDF, PDF], output_path=tmp_path)

    # Save a reference to the original pdf
    orig_pdf = API.pdf(pdf=PDF)

    pdf_alpha_s_bundle_name = PDF + '_pdfas'

    # Now try to read the PDF we just wrote to some non_std location
    with temp_lhapdf_path(tmp_path):
        pdf = lhapdf.mkPDFs(pdf_alpha_s_bundle_name)

        orig_pdf_members = orig_pdf.get_members()
        assert len(pdf) == orig_pdf_members + 2
        assert pdf[-1].type == pdf[-2].type == "central"

        # orig_pdf_members also includes replica0, therefore only +1
        assert pdf[-1].memberID == orig_pdf_members + 1
        assert pdf[-2].memberID == orig_pdf_members
