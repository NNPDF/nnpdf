#!/usr/bin/env python3
import pandas as pd
import numpy as np
from validphys.loader import FallbackLoader as Loader
from validphys.api import API
from collections import defaultdict
import sys
from nnpdf_data import legacy_to_new_map
import matplotlib.pyplot as plt
import warnings
import pathlib
from scipy.stats import norm
warnings.filterwarnings("ignore", category=RuntimeWarning)

def find_posterior(tcm_fitname, dest=pathlib.Path("/Users/s2850353/Documents/nnpdf_outputs/mc_determination")):
    fit = API.fit(fit=tcm_fitname)

    # find the name of point prescription mc_pp
    pps = fit.as_input()["theorycovmatconfig"]["point_prescriptions"]
    mc_pp_id, mc_pp = [[j,i] for j,i in enumerate(pps) if "mc" in i][0]

    # build validphys inputs
    common_dict = dict(
        dataset_inputs={"from_": "fit"},
        fit=fit.name,
        fits=[fit.name],
        use_cuts="fromfit",
        metadata_group="nnpdf31_process",
    )
    theoryids_dict = ({
            "point_prescription": mc_pp,
            "theoryid": {"from_": "theory"},
            "theory": {"from_": "fit"},
            "theorycovmatconfig": {"from_": "fit"},
        } | common_dict)
    
    # find theoryids
    theoryids = API.theoryids(**theoryids_dict)

    # goup theories by mc value. Potentially unstable.
    theory_plus = theoryids[2].id
    theory_mid = theoryids[0].id
    theory_min = theoryids[1].id

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

    # find priors
    prior_theorypreds_central = API.group_result_central_table_no_table(**inps_central)["theory_central"]
    prior_theorypreds_plus = API.group_result_central_table_no_table(**inps_plus)["theory_central"]
    prior_theorypreds_minus = API.group_result_central_table_no_table(**inps_minus)["theory_central"]

    # find the values of mc
    mc_plus = API.theory_info_table(theory_db_id=theory_plus).loc["mc"].iloc[0]
    mc_central = API.theory_info_table(theory_db_id=theory_mid).loc["mc"].iloc[0]
    mc_min = API.theory_info_table(theory_db_id=theory_min).loc["mc"].iloc[0]

    # make sure the mc shift in both directions is symmetric
    delta_mc_plus = mc_plus - mc_central
    delta_mc_min = mc_central - mc_min
    if abs(delta_mc_min - delta_mc_plus) > 1e-6:
        raise ValueError("mc shifts in both directions is not symmetric")
    else:
        mc_step_size = delta_mc_min

    # get the covmat scaling factor, it should be 1 in our case
    covmat_scaling_factor = fit.as_input().get("theorycovmatconfig",{}).get("rescale_alphas_covmat",1.0)

    # build the theory covmat
    beta_tilde = np.sqrt(1) * (mc_step_size / np.sqrt(2)) * np.array([1, -1])
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

    # read the fit covmat
    stored_mc_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_point_prescriptions{mc_pp_id}_theory_covmat_custom_per_prescription.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)


    storedcovmat_index = pd.MultiIndex.from_tuples(
        [(aa, bb, np.int64(cc)) for aa, bb, cc in stored_mc_covmat.index],
        names=["group", "dataset", "id"],
    )  # make sure theoryID is an integer, same as in S

    # format the covmat from fit
    stored_mc_covmat = pd.DataFrame(
        stored_mc_covmat.values, index=storedcovmat_index, columns=storedcovmat_index
    )
    new_names = {d[0]: legacy_to_new_map(d[0])[0] for d in stored_mc_covmat.index}
    stored_mc_covmat.rename(columns=new_names, index=new_names, level=1, inplace=True) # rename datasets using the legacy to new map
    stored_mc_covmat = stored_mc_covmat.reindex(S.index).T.reindex(S.index)

    # compare the covmats
    if not np.allclose(S, stored_mc_covmat):
        print("Reconstructed theory covmat, S, is not the same as the stored covmat!")

    # results per replica
    theorypreds_fit = API.group_result_table_no_table(**inps_central_fit).iloc[:, 2:]

    # experimental covmat
    exp_covmat = API.groups_covmat(
        use_t0=True,
        datacuts={"from_": "fit"},
        t0pdfset={"from_": "datacuts"},
        theoryid= {"from_": "theory"},
        theory={"from_": "fit"},
        **common_dict
    )

    # theory covmat
    total_th_covmat = pd.read_csv(
        fit.path / f"tables/datacuts_theory_theorycovmatconfig_theory_covmat_custom.csv",
        index_col=[0, 1, 2],
        header=[0, 1, 2],
        sep="\t|,",
        encoding="utf-8",
        engine="python",
    ).fillna(0)

    # reindex the theory covmat
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

    # mean prediction != prediction of mean PDF
    mean_prediction = theorypreds_fit.mean(axis=1)

    # calculate covariance matrix over replicas
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
    # pseudodata values per each replica
    dat_reps = dat_reps.reindex(S.index)

    # dat_central = API.group_result_central_table_no_table(**inps_central)["data_central"]
    # mean pseudodata (should be close to data)
    dat_central = dat_reps.mean(axis=1)

    # inverse of the full covmat
    invcov = pd.DataFrame(np.linalg.inv(exp_covmat + total_th_covmat),index=exp_covmat.index, columns=exp_covmat.index)
    invcov = invcov.reindex(S.index).T.reindex(S.index)

    # delta_T_tilde is Eq. 3.37 in https://arxiv.org/pdf/2105.05114
    # this is the shift in m_c
    delta_T_tilde = -S_hat @ invcov @ (mean_prediction - dat_central)

    # P_tilde is Eq. 3.38.
    # This is the variance. Need to include all the terms
    P_tilde = S_hat.T @ invcov @ X @ invcov @ S_hat + S_tilde - S_hat.T @ invcov @ S_hat

    # the result
    pred = mc_central + delta_T_tilde
    unc = np.sqrt(P_tilde)
    print(rf"Analytic prediction for $mc$: {pred:.5f} ± {unc:.5f}")

    # replica distribution
    tcm_vals = mc_central -S_hat @ invcov @ (theorypreds_fit.to_numpy() - dat_central.reindex(S.index).to_numpy().reshape(-1,1)) 

    # Adding artificial noise to approximate the replica distribution.
    noise = np.random.normal(loc=0.0, scale=1.0, size =tcm_vals.shape)
    stochastic_shift = np.sqrt(S_tilde-S_hat @ invcov @ S_hat.T) * noise
    tcm_vals+=stochastic_shift
    std_tcm = np.std(tcm_vals)
    mu_tcm, std_tcm = np.mean(tcm_vals), np.std(tcm_vals)
    print(f"Histogram mean and std: {mu_tcm:.5f} ± {std_tcm:.5f}")


    bins = np.linspace(1, 2, 30)
    # distribution of replicas
    plt.hist(tcm_vals, bins=bins, density=True, alpha=0.6, label='replica distribution', color="C0")

    x = np.linspace(bins[0], bins[-1], 1000)
    # fit Gaussian to replicas dist
    plt.plot(x, norm.pdf(x, pred, unc), '-', lw=2, color='orange', alpha=0.6, label='TCM')


    # Add labels and legend
    plt.xlabel(r'$m_\text{charm}$')
    plt.ylabel('normalized frequency')
    plt.title(f"{tcm_fitname}\nposterior mc={pred:.3f} ± {unc:.3f}")
    plt.tight_layout()
    plt.legend()
    filename=f"{tcm_fitname}_{pred}_pm_{unc}.pdf"
    plt.savefig(dest/filename)
    print("figure saved, all done.")

if __name__ == "__main__":
    tcm_fitname = sys.argv[1]
    if len(sys.argv) == 3: #script name and 2 args
        dest = pathlib.Path(sys.argv[2])
        find_posterior(tcm_fitname=tcm_fitname, dest=dest)
    else:
        find_posterior(tcm_fitname=tcm_fitname)