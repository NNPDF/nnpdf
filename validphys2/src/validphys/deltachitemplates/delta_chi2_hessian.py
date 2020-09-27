import numpy as np
from reportengine.checks import CheckError, make_argcheck

"""Check plot_delta_chi2_hessian is applied only to Hessian
set (symmhessian) converted from MC set.
"""


@make_argcheck
def check_pdf_is_symmhessian(pdf, **kwargs):
    etype = pdf.ErrorType
    if etype != "symmhessian":
        raise CheckError("Error: type of PDF %s must be 'symmhessian' and not %s" % (pdf, etype))


@check_pdf_is_symmhessian
def delta_chi2_hessian(pdf, experiments, groups_chi2):
    """
    Return delta_chi2 (computed as in plot_delta_chi2_hessian) relative to
    each eigenvector of the Hessian set.
    """
    experiments_chi2 = groups_chi2
    # store for each exp the chi2 from central value and chi2 from each error memebeer
    total_chis_exps = np.zeros(
        (len(experiments), 1 + len(experiments_chi2[0].replica_result.error_members()))
    )

    for i, ch in enumerate(experiments_chi2):
        th, central, _ = ch  # chi2: error member | central value | n data
        total_chis_exps[i] = [central, *th.error_members()]

    # sum over all exps to get chi2 total for cv and each error member
    total_chis = np.sum(total_chis_exps, axis=0)

    delta_chi2 = total_chis[1:] - total_chis[0]

    return delta_chi2
