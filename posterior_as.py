#!/usr/bin/env python3
import pandas as pd
import numpy as np
from validphys.loader import FallbackLoader as Loader
from validphys.api import API
from collections import defaultdict
import os
from nnpdf_data import legacy_to_new_map
import matplotlib.pyplot as plt
import warnings
from scipy.stats import norm
from pathlib import PosixPath
import sys
warnings.filterwarnings("ignore", category=RuntimeWarning)

def tcm_posterior(tcm_fitname, dest_dir):

    fit = API.fit(fit=tcm_fitname)

    # We have to know the name of the alphas point prescription (alphas_pp) to
    # extract the theoryids. We have to know alphas_pp_id to identify the .csv file
    # corresponding to the alphas covmat used in the fit
    pps = fit.as_input()["theorycovmatconfig"]["point_prescriptions"]
    alphas_pp_id, alphas_pp = [[j,i] for j,i in enumerate(pps) if "alphas" in i][0]




    common_dict = dict(
        dataset_inputs={"from_": "fit"},
        fit=fit.name,
        fits=[fit.name],
        use_cuts="fromfit",
        metadata_group="nnpdf31_process",
    )

    theoryids_dict = ({
            "point_prescription": alphas_pp,
            "theoryid": {"from_": "theory"},
            "theory": {"from_": "fit"},
            "theorycovmatconfig": {"from_": "fit"},
        } | common_dict)
    theoryids = API.theoryids(**theoryids_dict)
    theory_plus = theoryids[1].id
    theory_mid = theoryids[0].id
    theory_min = theoryids[2].id

    thcov_input_pdf = fit.as_input()["theorycovmatconfig"]["pdf"]




    # Inputs for central theory (used to construct the alphas covmat)
    inps_central = dict(theoryid=theory_mid, pdf=thcov_input_pdf, **common_dict)

    # Inputs for plus theory (used to construct the alphas covmat)
    inps_plus = dict(theoryid=theory_plus, pdf=thcov_input_pdf, **common_dict)

    # Inputs for minus theory prediction (used to construct the alphas covmat)
    inps_minus = dict(theoryid=theory_min, pdf=thcov_input_pdf, **common_dict)

    # inputs for the computation of the prediction of the fit with cov=C+S, where S
    # is computed using the inps_central, inps_plus, and inps_minus dictionaries
    inps_central_fit = dict(theoryid=theory_mid, pdf={"from_": "fit"}, **common_dict)


    # if fit.as_input()["resample_negative_pseudodata"] != False:
    #     print("The TCM assumes Gaussianity of the pseudodata, to ensure this set")
    #     print("resample_negative_pseudodata: False")
    #     print("in the n3fit runcard!")

    prior_theorypreds_central = API.group_result_central_table_no_table(**inps_central)["theory_central"]
    prior_theorypreds_plus = API.group_result_central_table_no_table(**inps_plus)["theory_central"]
    prior_theorypreds_minus = API.group_result_central_table_no_table(**inps_minus)["theory_central"]





    # Get the values of alphas...
    alphas_plus = API.theory_info_table(theory_db_id=theory_plus).loc["alphas"].iloc[0]
    alphas_central = API.theory_info_table(theory_db_id=theory_mid).loc["alphas"].iloc[0]
    alphas_min = API.theory_info_table(theory_db_id=theory_min).loc["alphas"].iloc[0]

    # ... and make sure the alphas shift in both directions is symmetric
    delta_alphas_plus = alphas_plus - alphas_central
    delta_alphas_min = alphas_central - alphas_min
    if abs(delta_alphas_min - delta_alphas_plus) > 1e-6:
        raise ValueError("alphas shifts in both directions is not symmetric")
    else:
        alphas_step_size = delta_alphas_min





    # At some point we scaled the covmat to account for higher order derivatives or
    # to test depencence of the prior. It is not used in the final result
    covmat_scaling_factor = fit.as_input().get("theorycovmatconfig",{}).get("rescale_alphas_covmat",1.0)

    beta_tilde = np.sqrt(1) * (alphas_step_size / np.sqrt(2)) * np.array([1, -1])
    S_tilde = beta_tilde @ beta_tilde


    delta_plus = (np.sqrt(covmat_scaling_factor) / np.sqrt(2)) * (
        prior_theorypreds_plus - prior_theorypreds_central
    )
    delta_minus = (np.sqrt(covmat_scaling_factor) / np.sqrt(2)) * (
        prior_theorypreds_minus - prior_theorypreds_central
    )

    beta = [delta_plus, delta_minus]
    S_hat = pd.Series(beta_tilde @ beta, index=delta_minus.index)

    S = np.outer(delta_plus, delta_plus) + np.outer(delta_minus, delta_minus)
    S = pd.DataFrame(S, index=delta_minus.index, columns=delta_minus.index)


    stored_alphas_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_point_prescriptions{alphas_pp_id}_theory_covmat_custom_per_prescription.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)
    storedcovmat_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in stored_alphas_covmat.index],
        names=["group", "dataset", "id"],
    )  # make sure theoryID is an integer, same as in S
    stored_alphas_covmat = pd.DataFrame(
        stored_alphas_covmat.values, index=storedcovmat_index, columns=storedcovmat_index
    )
    new_names = {d[0]: legacy_to_new_map(d[0])[0] for d in stored_alphas_covmat.index}
    stored_alphas_covmat.rename(columns=new_names, index=new_names, level=1, inplace=True) # rename datasets using the legacy to new map
    stored_alphas_covmat = stored_alphas_covmat.reindex(S.index).T.reindex(S.index)
    # stored_alphas_covmat = stored_alphas_covmat.reindex(S.index).T.reindex(S.index)

    if not np.allclose(S, stored_alphas_covmat):
        print("Reconstructed theory covmat, S, is not the same as the stored covmat!")


    theorypreds_fit = API.group_result_table_no_table(**inps_central_fit).iloc[:, 2:]


    exp_covmat = API.groups_covmat(
        use_t0=True,
        datacuts={"from_": "fit"},
        t0pdfset={"from_": "datacuts"},
        theoryid= {"from_": "theory"},
        theory={"from_": "fit"},
        **common_dict
    )


    total_th_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_theory_covmat_custom.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)
    new_names = {d[0]: legacy_to_new_map(d[0])[0] for d in total_th_covmat.index}
    total_th_covmat.rename(columns=new_names, index=new_names, level=1, inplace=True) # rename datasets using the legacy to new map
    total_th_covmat_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in total_th_covmat.index],
        names=["group", "dataset", "id"],
    ) # make sure the index is an int, just as it is in S
    total_th_covmat = pd.DataFrame(
        total_th_covmat.values, index=total_th_covmat_index, columns=total_th_covmat_index
    )
    total_th_covmat = total_th_covmat.reindex(S.index).T.reindex(S.index)


    # Note that mean_prediction is different from the prediction of the mean PDF (i.e. replica0)
    mean_prediction = theorypreds_fit.mean(axis=1)

    X = np.zeros_like(S.values)
    for i in range(theorypreds_fit.shape[1]):
        X += np.outer(
            (theorypreds_fit.iloc[:, i] - mean_prediction),
            (theorypreds_fit.iloc[:, i] - mean_prediction),
        )
    X *= 1 / theorypreds_fit.shape[1]
    X = pd.DataFrame(X, index=theorypreds_fit.index, columns=theorypreds_fit.index)


    # In the computation we use <D>_rep and not the central value of the data D_exp, though if
    # resample_negative_pseudodata: false
    # is set in the n3fit runcard, D_exp and <D>_rep should be the same as N_rep -> inf.
    pseudodata = API.read_pdf_pseudodata(**common_dict)
    dat_reps = pd.concat(
        [i.pseudodata for i in pseudodata], axis=1
    )
    dat_reps = dat_reps.reindex(S.index)


    # dat_central = API.group_result_central_table_no_table(**inps_central)["data_central"]
    dat_central = dat_reps.mean(axis=1)


    invcov = pd.DataFrame(np.linalg.inv(exp_covmat + total_th_covmat),index=exp_covmat.index, columns=exp_covmat.index)
    invcov = invcov.reindex(S.index).T.reindex(S.index)


    # delta_T_tilde is Eq. 3.37 in https://arxiv.org/pdf/2105.05114
    delta_T_tilde = -S_hat @ invcov @ (mean_prediction - dat_central)

    # P_tilde is Eq. 3.38.
    #
    # Note that not all terms of the equation in the paper are here, in particular
    # X_tile and X_hat vanish. This is because they measure the covariance of
    # T_tilde over PDF replicas, but for us T_tilde is alphas. The prediciton of
    # alphas does not depend on the PDF, and as such T_tilde^(r) == T_tilde^(0)
    P_tilde = S_hat.T @ invcov @ X @ invcov @ S_hat + S_tilde - S_hat.T @ invcov @ S_hat

    pred = alphas_central + delta_T_tilde
    unc = np.sqrt(P_tilde)
    print(rf"Prediction for $\alpha_s$: {pred:.5f} ± {unc:.5f}")

    # =============== REPLICA DISTRIBUTION ====================================
    tcm_vals = alphas_central -S_hat @ invcov @ (theorypreds_fit.to_numpy() - dat_central.reindex(S.index).to_numpy().reshape(-1,1))

    bins = np.linspace(0.117, 0.122, 31)
    plt.hist(tcm_vals, bins=bins, density=True, alpha=0.6, label='TCM', color="C0")

    x = np.linspace(bins[0], bins[-1], 1000)
    mu_tcm, std_tcm = np.mean(tcm_vals), np.std(tcm_vals)
    plt.plot(x, norm.pdf(x, mu_tcm, std_tcm), 'b-', lw=2, color="C0", alpha=0.6)


    # Add labels and legend
    plt.xlabel(r'$\alpha_s(m_Z)$')
    plt.ylabel('normalized frequency')
    # title = "with positivity" if positivity else "without positivity"
    plt.title(r"$\mathrm{NNLO}_\mathrm{QCD}$")
    plt.tight_layout()
    plt.legend()
    plt.show()
    figname=f"{tcm_fitname}_{pred:.5f}_pm_{unc:.5f}".replace(".", "p")+".pdf"
    dest_dir=PosixPath(dest_dir)
    plt.savefig(dest_dir / figname)

if __name__ == "__main__":
    tcm_fitname = sys.argv[1]
    dest_dir = sys.argv[2]
    tcm_posterior(tcm_fitname=tcm_fitname, dest_dir=dest_dir)
