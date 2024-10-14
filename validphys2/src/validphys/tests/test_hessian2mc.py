import numpy as np
import pandas as pd
from unittest import mock
from validphys.hessian2mc import write_mc_watt_thorne_replicas, write_hessian_to_mc_watt_thorne
import pathlib


@mock.patch("validphys.hessian2mc.write_replica")
@mock.patch("validphys.hessian2mc.log.info")
def test_write_mc_watt_thorne_replicas(mock_log_info, mock_write_replica):
    # Set up mock inputs
    n_eig = 5
    n_all_eig = 1 + 2 * n_eig  # includes central member, so 1 + 2*n_eig
    num_members = 50  # new MC members
    num_nx_nq = 100  # xgrid and Q points

    replicas_df = pd.DataFrame(
        np.random.standard_normal(size=(num_nx_nq, n_all_eig)), columns=range(1, n_all_eig + 1)
    )

    Rjk_std_normal = np.random.standard_normal(size=(num_members, n_eig))

    # import IPython; IPython.embed()
    mc_pdf_path = mock.Mock()  # Mock the path

    # Call the function being tested
    write_mc_watt_thorne_replicas(Rjk_std_normal, replicas_df, mc_pdf_path)

    # Verify that write_replica was called the correct number of times
    assert (
        mock_write_replica.call_count == num_members + 1
    )  # for Rjk members + 1 for the central replica

    # Check if the log message was correct for the central replica
    mock_log_info.assert_any_call(f"Writing central replica to {mc_pdf_path}")


@mock.patch("validphys.hessian2mc.write_mc_watt_thorne_replicas")
@mock.patch("validphys.hessian2mc.load_all_replicas")
@mock.patch("validphys.hessian2mc.rep_matrix")
@mock.patch("validphys.hessian2mc.write_new_lhapdf_info_file_from_previous_pdf")
@mock.patch("validphys.hessian2mc.os.makedirs")
@mock.patch("validphys.hessian2mc.lhapdf.paths")
@mock.patch("validphys.hessian2mc.PDF")
def test_write_hessian_to_mc_watt_thorne(
    mock_pdf,
    mock_lhapdf_paths,
    mock_makedirs,
    mock_write_info_file,
    mock_rep_matrix,
    mock_load_all_replicas,
    mock_write_mc_replicas,
):
    # Set up mock inputs
    msht_like_hessian_pdf = "MSHT20"
    mc_pdf_name = "MSHT20_MC"
    num_members = 100

    mock_load_all_replicas.return_value = (None, None)

    mock_pdf_instance = mock.Mock()
    mock_pdf.return_value = mock_pdf_instance

    mock_lhapdf_paths.return_value = [pathlib.Path("/path/to/lhapdf")]

    mock_rep_matrix.return_value = np.random.randn(5, 7)  # Mocked replica matrix

    # Call the function being tested
    write_hessian_to_mc_watt_thorne(msht_like_hessian_pdf, mc_pdf_name, num_members)

    # Verify that the necessary directories were created
    mc_pdf_path = pathlib.Path("/path/to/lhapdf") / mc_pdf_name
    mock_makedirs.assert_called_once_with(mc_pdf_path)

    # Verify that the LHAPDF info file was written
    mock_write_info_file.assert_called_once_with(
        path_old_pdfset=pathlib.Path("/path/to/lhapdf") / msht_like_hessian_pdf,
        name_old_pdfset=msht_like_hessian_pdf,
        path_new_pdfset=mc_pdf_path,
        name_new_pdfset=mc_pdf_name,
        num_members=num_members,
        description_set=f"MC representation of {msht_like_hessian_pdf}",
        errortype="replicas",
    )

    # Verify that the replicas were written
    mock_write_mc_replicas.assert_called_once()
