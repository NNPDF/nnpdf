#!/usr/bin/env python3
import numpy as np
import pandas as pd
from validphys.api import API
from nnpdf_data import legacy_to_new_map
import sys
import pathlib

import matplotlib.pyplot as plt
from matplotlib import patches, transforms
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Ellipse
import scipy


def theory_cov_method(fitname):
    fit = API.fit(fit=fitname)


    # find point prescription names
    pps = fit.as_input()["theorycovmatconfig"]["point_prescriptions"]
    for i, pp in enumerate(pps):
            if "mcharm" in pp:
                mcharm_pp_id = i
                mcharm_pp = pp
            elif "alphas" in pp:
                alphas_pp_id = i
                alphas_pp = pp
            else:
                continue


    common_dict = dict(
            dataset_inputs={"from_": "fit"},
            fit=fit.name,
            fits=[fit.name],
            use_cuts="fromfit",
            metadata_group="nnpdf31_process",
        )


    theoryids_dict_mcharm = ({
                                "point_prescription": mcharm_pp,
                                "theoryid": {"from_": "theory"},
                                "theory": {"from_": "fit"},
                                "theorycovmatconfig": {"from_": "fit"},
                            } | common_dict)


    theoryids_dict_alphas = ({
                                    "point_prescription": alphas_pp,
                                    "theoryid": {"from_": "theory"},
                                    "theory": {"from_": "fit"},
                                    "theorycovmatconfig": {"from_": "fit"},
                                } | common_dict)


    # extract theory ids for mcharm and alphas including their variations
    theoryids_mcharm = API.theoryids(**theoryids_dict_mcharm)
    theoryids_alphas = API.theoryids(**theoryids_dict_alphas)


    theory_central = theoryids_alphas[0].id #=theoryids_mcharm[0].id
    theory_plus = [theoryids_mcharm[2].id, theoryids_alphas[1].id]
    theory_min = [theoryids_mcharm[1].id, theoryids_alphas[2].id]


    thcov_input_pdf = fit.as_input()["theorycovmatconfig"]["pdf"]


    # Inputs for central theory (used to construct the mcharm, alphas covmat).
    inps_central = dict(theoryid=theory_central, pdf=thcov_input_pdf, **common_dict)

    # Inputs for plus theory (used to construct the mcharm, alphas covmat)
    inps_plus = [dict(theoryid=theory_id_plus, pdf=thcov_input_pdf, **common_dict) for theory_id_plus in theory_plus]

    # Inputs for minus theory prediction (used to construct the mcharm, alphas covmat)
    inps_min = [dict(theoryid=theory_id_min, pdf=thcov_input_pdf, **common_dict) for theory_id_min in theory_min]



    # inputs for the computation of the prediction of the fit with cov=C+S, where S
    # is computed using the inps_central, inps_plus, and inps_minus dictionaries
    inps_central_fit = dict(theoryid=theory_central, pdf={"from_": "fit"}, **common_dict)



    # if fit.as_input()["resample_negative_pseudodata"] != False:
    #     print("The TCM assumes Gaussianity of the pseudodata, to ensure this set")
    #     print("resample_negative_pseudodata: False")
    #     print("in the n3fit runcard!")


    # Get the prior theory predictions for the central values
    prior_theorypreds_central = API.group_result_central_table_no_table(**inps_central)[
        "theory_central"].to_frame()  # shape (n_dat, 1)

    # Get the prior theory predictions for the plus and minus variations
    prior_theorypreds_plus = pd.concat(
        [API.group_result_central_table_no_table(**inp_plus)["theory_central"] for inp_plus in inps_plus],
        axis=1)  # shape (n_dat, n_par)
    prior_theorypreds_minus = pd.concat(
        [API.group_result_central_table_no_table(**inp_min)["theory_central"] for inp_min in inps_min],
        axis=1)  # shape (n_dat, n_par)



    # Get the values of mcharm
    mcharm_plus = API.theory_info_table(theory_db_id=theory_plus[0]).loc["mc"].iloc[0]
    mcharm_central = API.theory_info_table(theory_db_id=theory_central).loc["mc"].iloc[0]
    mcharm_min = API.theory_info_table(theory_db_id=theory_min[0]).loc["mc"].iloc[0]

    # and alphas
    alphas_plus = API.theory_info_table(theory_db_id=theory_plus[1]).loc["alphas"].iloc[0]
    alphas_central = API.theory_info_table(theory_db_id=theory_central).loc["alphas"].iloc[0]
    alphas_min = API.theory_info_table(theory_db_id=theory_min[1]).loc["alphas"].iloc[0]



    # and make sure the shift in both directions is symmetric
    delta_plus = np.array([mcharm_plus - mcharm_central, alphas_plus - alphas_central])
    delta_min = np.array([mcharm_central - mcharm_min, alphas_central - alphas_min])
    if np.any(abs(delta_min - delta_plus) > 1e-6):
        raise ValueError("mcharm shifts in both directions is not symmetric")
    else:
        step_size = np.array(delta_min).reshape(-1, 1)


    # At some point we scaled the covmat to account for higher order derivatives or
    # to test depencence of the prior. It is not used in the final result
    covmat_scaling_factor = 1  # fit.as_input().get("theorycovmatconfig",{}).get("rescale_alphas_covmat",1.0)
    lambda_index = ["mcharm", "alphas"]


    # Compute theory covmat S_tilde on the genuine predictions (as, mt)
    beta_tilde = np.sqrt(covmat_scaling_factor) * (step_size / np.sqrt(2)) * np.array([[1, -1, 0, 0], [0, 0, 1, -1]])
    S_tilde = pd.DataFrame(beta_tilde @ beta_tilde.T, index=lambda_index, columns=lambda_index)

    # Overriding a variable!
    delta_plus = (np.sqrt(covmat_scaling_factor) / np.sqrt(2)) * (
            prior_theorypreds_plus - prior_theorypreds_central
    )
    delta_minus = (np.sqrt(covmat_scaling_factor) / np.sqrt(2)) * (
            prior_theorypreds_minus - prior_theorypreds_central
    )


    # Compute the theory cross covmat between the fitted predictions and the genuine predictions
    beta = np.array([delta_plus.iloc[:, 0].values, delta_minus.iloc[:, 0].values, delta_plus.iloc[:, 1].values,
                        delta_minus.iloc[:, 1].values]).T  # shape (n_dat, 2 * n_par)
    S_hat = pd.DataFrame(beta_tilde @ beta.T, columns=delta_minus.index, index=lambda_index)  # shape (n_par, n_dat)

    # Compute the theory covmat on the theory predictions
    S = pd.DataFrame(beta @ beta.T, index=delta_minus.index, columns=delta_minus.index)


    # experimental covmat
    exp_covmat = API.groups_covmat(
        use_t0=True,
        datacuts={"from_": "fit"},
        t0pdfset={"from_": "datacuts"},
        theoryid={"from_": "theory"},
        theory={"from_": "fit"},
        **common_dict
    )


    stored_alphas_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_point_prescriptions{alphas_pp_id}_theory_covmat_custom_per_prescription.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)


    stored_mcharm_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_point_prescriptions{mcharm_pp_id}_theory_covmat_custom_per_prescription.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)


    stored_theory_covmat_all = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_theory_covmat_custom.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)


    # convert 3rd index to int
    stored_alphas_covmat_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in stored_alphas_covmat.index],
        names=["group", "dataset", "id"],
    )

    stored_mcharm_covmat_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in stored_mcharm_covmat.index],
        names=["group", "dataset", "id"],
    )

    stored_theory_covmat_all_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in stored_theory_covmat_all.index],
        names=["group", "dataset", "id"],
    )


    # make sure theoryID is an integer, same as in S
    stored_alphas_covmat = pd.DataFrame(
        stored_alphas_covmat.values, index=stored_alphas_covmat_index, columns=stored_alphas_covmat_index
    )

    stored_mcharm_covmat = pd.DataFrame(
        stored_mcharm_covmat.values, index=stored_mcharm_covmat_index, columns=stored_mcharm_covmat_index
    )

    stored_theory_covmat_all = pd.DataFrame(stored_theory_covmat_all.values, index=stored_theory_covmat_all_index,
                                            columns=stored_theory_covmat_all_index)

    new_names = {d[1]: legacy_to_new_map(d[1])[0] for d in stored_alphas_covmat.index}



    # rename datasets using the legacy to new map
    stored_alphas_covmat.rename(columns=new_names, index=new_names, level=1, inplace=True)
    stored_mcharm_covmat.rename(columns=new_names, index=new_names, level=1, inplace=True)
    stored_theory_covmat_all.rename(columns=new_names, index=new_names, level=1, inplace=True)

    stored_alphas_covmat = stored_alphas_covmat.reindex(S.index).T.reindex(S.index)
    stored_mcharm_covmat = stored_mcharm_covmat.reindex(S.index).T.reindex(S.index)
    stored_theory_covmat_all = stored_theory_covmat_all.reindex(S.index).T.reindex(S.index)

    alphas_mcharm_covmat = stored_alphas_covmat + stored_mcharm_covmat
    if not np.allclose(S, alphas_mcharm_covmat):
        print("Reconstructed theory covmat, S, is not the same as the stored covmat!")
    if np.allclose(S.to_numpy(), alphas_mcharm_covmat.to_numpy()):
        print("values are close.")
    # aren't we forgetting about MHOUs? How do they work?
    data_theory_results = API.group_result_table_no_table(**inps_central_fit)
    theorypreds_fit = data_theory_results.iloc[:, 2:]



    # Note that mean_prediction is different from the prediction of the mean PDF (i.e. replica0)
    T0 = theorypreds_fit.mean(axis=1)
    X = np.cov(theorypreds_fit)



    # X, _ = remove_negative_eigenmodes(X, tol=0)
    X = pd.DataFrame(X, index=theorypreds_fit.index, columns=theorypreds_fit.index)


    # remove negative eigenmodes from X

    # In the computation we use <D>_rep and not the central value of the data D_exp, though if
    # resample_negative_pseudodata: false
    # is set in the n3fit runcard, D_exp and <D>_rep should be the same as N_rep -> inf.
    pseudodata = API.read_pdf_pseudodata(**common_dict)
    dat_reps = pd.concat(
        [i.pseudodata.reindex(prior_theorypreds_central.index) for i in pseudodata], axis=1
    )

    invcov = pd.DataFrame(np.linalg.inv(exp_covmat + stored_theory_covmat_all), index=exp_covmat.index,
                            columns=exp_covmat.index)


    # delta_T_tilde = -S_hat @ invcov @ (T0 - dat_reps.mean(axis=1))
    delta_T_tilde = -S_hat @ invcov @ (T0 - data_theory_results["data_central"])


    # P_tilde is Eq. 3.38.
    #
    # Note that not all terms of the equation in the paper are here, in particular
    # X_tile and X_hat vanish. This is because they measure the covariance of
    # T_tilde over PDF replicas, but for us T_tilde is alphas. The prediciton of
    # alphas does not depend on the PDF, and as such T_tilde^(r) == T_tilde^(0)

    P_tilde = S_hat @ invcov @ X @ invcov @ S_hat.T + S_tilde - S_hat @ invcov @ S_hat.T
    central_theory = np.array([mcharm_central, alphas_central])
    pred = central_theory + delta_T_tilde

    str_result = print_results(pred, P_tilde)
    return pred, P_tilde, str_result



def remove_negative_eigenmodes(matrix, tol=1e-12):
    """
    Removes negative eigenmodes from a matrix by projecting onto the
    subspace of non-negative eigenvalues.

    Parameters
    ----------
    matrix : (n, n) array_like
        Symmetric (or Hermitian) matrix.
    tol : float
        Tolerance below which eigenvalues are considered zero.

    Returns
    -------
    matrix_pos : (n, n) ndarray
        Matrix reconstructed with only non-negative eigenmodes.
    eigvals : ndarray
        Eigenvalues after truncation (negative ones set to 0).
    """
    # Ensure the matrix is numpy array
    A = np.array(matrix, dtype=float)

    # Compute eigen-decomposition
    eigvals, eigvecs = np.linalg.eigh(A)
    print(f"number of negative evals in X is:  {sum(eigvals<tol)} out of {len(eigvals)}")
    # Zero out negative eigenvalues (or those below tolerance)
    eigvals_clipped = np.where(eigvals > tol, eigvals, 0.0)

    # Reconstruct the positive semidefinite matrix
    A_pos = eigvecs @ np.diag(eigvals_clipped) @ eigvecs.T

    return A_pos, eigvals_clipped


def print_results(central_value, cov):
    central_mcharm = central_value.loc["mcharm"]
    central_alphas = central_value.loc["alphas"]
    sigma_mcharm = np.sqrt(cov["mcharm"]["mcharm"])
    sigma_alphas = np.sqrt(cov["alphas"]["alphas"])
    rho = cov["mcharm"]["alphas"] / (sigma_mcharm * sigma_alphas)
    str_result = f"mcharm (68%): {central_mcharm:.3f} +/- {sigma_mcharm:.3f} \n alphas (68%): {central_alphas:.6f} +/- {sigma_alphas:.6f} \n rho: {rho:.3f}"
    print(str_result)
    return str_result

def confidence_ellipse(ax, cov, mean, facecolor=None, confidence_level=95, **kwargs):
    handles = []
    # find eigs of covmat
    eig_val, eig_vec = np.linalg.eigh(cov)
    order = np.argsort(eig_val)[::-1]
    eig_val, eig_vec = eig_val[order], eig_vec[:, order]

    # find s
    chi2_qnt = scipy.stats.chi2.ppf(confidence_level / 100.0, 2)
    width, height = 2 * np.sqrt(chi2_qnt * eig_val)
    angle = np.degrees(np.arctan2(*eig_vec[:, 0][::-1]))

    mean = np.asarray(mean, dtype=float).ravel()

    ellipse = Ellipse(
        xy=mean, width=width, height=height, angle=angle,
        facecolor=facecolor, edgecolor=kwargs.get("edgecolor", "black"), linewidth=kwargs.get("linewidth", 1.5),
        linestyle=kwargs.get("linestyle", "-"),
    )

    if facecolor is not None:
        ellipse.set_facecolor((*plt.cm.colors.to_rgba(facecolor)[:3], 0.15))
    else:
        ellipse.set_facecolor("none")
    ellipse.set_edgecolor((*plt.cm.colors.to_rgba(kwargs.get("edgecolor", "black"))[:3], 1.0))
    handles.append(ellipse)
    ax.add_patch(ellipse)
    if kwargs.get("marker", True):
        ax.scatter(mean[0], mean[1], marker="x", color="black")
    ax.grid(True, which='major', linestyle='-', linewidth=0.8)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))  # number of minor intervals per major tick
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which='minor', linestyle=':', linewidth=0.6, alpha=0.7)
    return width, height, handles
