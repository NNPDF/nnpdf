"""
    Test for the mc2hessian module
"""
from validphys.api import API
from validphys.tests.conftest import temp_lhapdf_path

NEIG = 5

def test_mc2hessian(data_config, tmp):
    """Tests that the generated hessian PDF is indeed marked as such
    and that the metadata is not obviously broken
    """
    config = dict(data_config)
    pdf_hessian_name = f"{config['pdf']}_hessian_{NEIG}"

    config["Neig"] = NEIG
    config["output_path"] = tmp
    API.mc2hessian(**config)

    # Save a reference to the original pdf
    orig_pdf_info = API.pdf(pdf=config["pdf"]).info

    # Now try to read the PDF we just wrote to some non_std location
    with temp_lhapdf_path(tmp.as_posix()):
        pdf = API.pdf(pdf=pdf_hessian_name)

        # Check the correct changed metadata
        assert pdf.error_type == "symmhessian"
        assert len(pdf) == NEIG + 1

        # Now check that all info is the same as the original
        # except for ErrorType that _must_ be different
        skip_this = ["SetDesc", "NumMembers"]
        for key, item in orig_pdf_info.items():
            if key in skip_this:
                continue
            new_item = pdf.info[key]
            if key == "ErrorType":
                assert item != new_item
            else:
                assert item == new_item
