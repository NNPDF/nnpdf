"""
validphys.hessian2mc.py

This module contains the functions that can be used to convert Hessian sets
like MSHT20 and CT18 to Monte Carlo sets.
The functions implemented here follow equations (4.3), for MSHT20, and (4.4), for CT18,
of the paper arXiv:2203.05506
"""

import pathlib
import lhapdf
import os
import logging
import numpy as np

from validphys.lhio import load_all_replicas, rep_matrix, write_replica
from validphys.checks import check_pdf_is_hessian

log = logging.getLogger(__name__)


def write_new_lhapdf_info_file_from_previous_pdf(
    path_old_pdfset,
    name_old_pdfset,
    path_new_pdfset,
    name_new_pdfset,
    num_members,
    description_set="MC representation of hessian PDF set",
    errortype="replicas",
):
    """
    Writes a new LHAPDF set info file based on an existing set.
    """

    # write LHAPDF info file for a new pdf set
    with open(path_old_pdfset / f"{name_old_pdfset}.info", "r") as in_stream, open(
        path_new_pdfset / f"{name_new_pdfset}.info", "w"
    ) as out_stream:
        for l in in_stream.readlines():
            if l.find("SetDesc:") >= 0:
                out_stream.write(f'SetDesc: f"{description_set}"\n')
            elif l.find("NumMembers:") >= 0:
                out_stream.write(f"NumMembers: {num_members}\n")
            elif l.find("ErrorType:") >= 0:
                out_stream.write(f"ErrorType: {errortype}\n")
            elif l.find("ErrorConfLevel") >= 0:
                # remove ErrorConfLevel line
                pass
            else:
                out_stream.write(l)
    log.info(f"Info file written to {path_new_pdfset / f'{name_new_pdfset}.info'}")


def write_mc_watt_thorne_replicas(Rjk_std_normal, replicas_df, mc_pdf_path):
    """
    Writes the Monte Carlo representation of a PDF set that is in Hessian form
    using the Watt-Thorne (MSHT20) prescription described in Eq. 4.3 of arXiv:2203.05506.

    Parameters
    ----------
    Rjk_std_normal: np.ndarray
        Array of shape (num_members, n_eig) containing random standard normal numbers.

    replicas_df: pd.DataFrame
        DataFrame containing replicas of the hessian set at all scales.

    mc_pdf_path: pathlib.Path
        Path to the new Monte Carlo PDF set.
    """

    for i, rnd_std_norm_vec in enumerate(Rjk_std_normal):

        # Odd eigenvectors: negative direction, even eigenvectors: positive direction
        df_odd = replicas_df.loc[:, 2::2]
        df_even = replicas_df.loc[:, 3::2]
        new_column_names = range(1, len(df_even.columns) + 1)

        df_even.columns = new_column_names
        df_odd.columns = new_column_names

        central_member, hess_diff_cov = replicas_df.loc[:, [1]], df_even - df_odd

        # Eq. 4.3 of arXiv:2203.05506
        mc_replica = central_member.dot([1]) + 0.5 * hess_diff_cov.dot(rnd_std_norm_vec)

        wm_headers = f"PdfType: replica\nFormat: lhagrid1\nFromMCReplica: {i}\n"
        log.info(f"Writing replica {i + 1} to {mc_pdf_path}")
        write_replica(i + 1, mc_pdf_path, wm_headers.encode("UTF-8"), mc_replica)

    # Write central replica from hessian set to mc set
    wm_headers = f"PdfType: replica\nFormat: lhagrid1\nFromMCReplica: {i}\n"
    log.info(f"Writing central replica to {mc_pdf_path}")
    write_replica(0, mc_pdf_path, wm_headers.encode("UTF-8"), central_member)


@check_pdf_is_hessian
def write_hessian_to_mc_watt_thorne(pdf, mc_pdf_name, num_members, watt_thorne_rnd_seed=1):
    """
    Writes the Monte Carlo representation of a PDF set that is in Hessian form
    using the Watt-Thorne (MSHT20) prescription described in Eq. 4.3 of arXiv:2203.05506.

    Parameters
    ----------
    pdf: validphys.core.PDF
        The Hessian PDF set that is to be converted to Monte Carlo.

    mc_pdf_name: str
        The name of the new Monte Carlo PDF set.

    """
    hessian_set = pdf

    lhapdf_path = pathlib.Path(lhapdf.paths()[-1])

    # path to hessian lhapdf set
    hessian_pdf_path = lhapdf_path / str(hessian_set)

    # path to new wmin pdf set
    mc_pdf_path = lhapdf_path / mc_pdf_name

    # create new wmin pdf set folder in lhapdf path if it does not exist
    if not mc_pdf_path.exists():
        os.makedirs(mc_pdf_path)

    # write LHAPDF info file for new wmin pdf set
    write_new_lhapdf_info_file_from_previous_pdf(
        path_old_pdfset=hessian_pdf_path,
        name_old_pdfset=str(hessian_set),
        path_new_pdfset=mc_pdf_path,
        name_new_pdfset=mc_pdf_name,
        num_members=num_members,
        description_set=f"MC representation of {str(hessian_set)}",
        errortype="replicas",
    )

    # load replicas from basis set at all scales
    _, grids = load_all_replicas(hessian_set)
    replicas_df = rep_matrix(grids)

    # Eg for MSHT20, mem=0 => central value; mem=1-64 => 32 eigenvector sets (+/- directions)
    n_eig = int((replicas_df.shape[1] - 1) / 2)

    np.random.seed(watt_thorne_rnd_seed)
    Rjk_std_normal = np.random.standard_normal(size=(num_members, n_eig))

    # write replicas to new wmin pdf set
    write_mc_watt_thorne_replicas(Rjk_std_normal, replicas_df, mc_pdf_path)
