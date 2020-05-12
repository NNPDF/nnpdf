"""
closuretest/checks.py

Module containing checks specfic to the closure tests
"""
from reportengine.checks import make_argcheck, CheckError


@make_argcheck
def check_use_fitcommondata(use_fitcommondata):
    """Base check that `use_fitcommondata` is being used, check should be used
    with all actions which require comparison to fitcommondata
    """
    if not use_fitcommondata:
        raise CheckError(
            "use_fitcommondata must be set to True for closure test estimators"
        )


@make_argcheck
def check_fit_isclosure(fit):
    """Check the input fit is a closure test"""
    if not fit.as_input()["closuretest"]["fakedata"]:
        raise CheckError(f"Specified fit: {fit}, is not a closure test")


@make_argcheck
def check_fits_underlying_law_match(fits):
    """Check that the fits being compared have the same underlying law"""
    # check that set of underlying laws is single item
    laws = {fit.as_input()["closuretest"]["fakepdf"] for fit in fits}
    if len(laws) != 1:
        raise CheckError(
            f"Closure test fits were fitting pseudo data generated from different input PDFs: {laws}"
        )


@make_argcheck
def check_fits_same_filterseed(fits):
    """Input fits should have the same filter seed if they are being compared"""
    seeds = {fit.as_input()["closuretest"]["filterseed"] for fit in fits}
    if len(seeds) != 1:
        raise CheckError(
            f"Closure test fits were fitting pseudo data generated with different level 1 noise: {seeds}"
        )


@make_argcheck
def check_fits_areclosures(fits):
    """Check all fits are closures"""
    for fit in fits:
        if not fit.as_input()["closuretest"]["fakedata"]:
            raise CheckError(f"Specified fit: {fit}, is not a closure test")


@make_argcheck
def check_t0_used(use_t0):
    if not use_t0:
        raise CheckError("use_t0 must be true")


@make_argcheck
def check_t0pdfset_matches_law(t0pdfset, fit):
    t0_from_fit = fit.as_input()["closuretest"]["fakepdf"]
    if not str(t0pdfset) == t0_from_fit:
        raise CheckError(
            f"Underlying pdf: {t0_from_fit}, does not match t0pdfset: {t0pdfset}"
        )


@make_argcheck
def check_at_least_10_fits(fits):
    if len(fits) < 10:
        raise CheckError(
            "Multiclosure actions testing finite sampling effects require at least 10 fits"
        )


@make_argcheck
def check_multifit_replicas(fits_pdf, _internal_n_reps):
    """Checks that all the fit pdfs have the same number of replicas N_rep and
    that N_rep is at least 10"""
    # we take off 1 here because we don't want to include replica 0
    if _internal_n_reps is not None:
        return {"_internal_n_reps": _internal_n_reps}

    n_reps = set([len(pdf) - 1 for pdf in fits_pdf])
    if len(n_reps) != 1:
        raise CheckError(
            "all fits for multiclosure actions should have same number of replicas"
        )
    n_reps = n_reps.pop()
    if n_reps < 10:
        raise CheckError(
            "Multiclosure actions testing finite sampling effects require fits "
            "to have at least 10 replicas"
        )
    return {"_internal_n_reps": n_reps}


@make_argcheck
def check_fits_different_filterseed(fits):
    """Input fits should have the same filter seed if they are being compared"""
    seed_list = [fit.as_input()["closuretest"]["filterseed"] for fit in fits]
    if len(seed_list) > len(set(seed_list)):
        raise CheckError(
            f"Multiclosure actions require that fits had different level 1 noise"
        )
