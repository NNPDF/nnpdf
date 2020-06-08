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
def check_multifit_replicas(fits_pdf, _internal_max_reps, _internal_min_reps):
    """Checks that all the fit pdfs have the same number of replicas N_rep and
    that N_rep is at least 20 (by default).

    This check also has the secondary
    effect of filling in the namespace keys _internal_max_reps and _internal_min_reps
    which can be used to override the number of replicas used at the level of
    the runcard, but by default get filled in as the number of replicas in each
    fit and 20 respectively

    """
    # we take off 1 here because we don't want to include replica 0
    n_reps = {len(pdf) - 1 for pdf in fits_pdf}
    if len(n_reps) != 1:
        raise CheckError(
            "all fits for multiclosure actions should have same number of replicas"
        )
    n_reps = n_reps.pop()
    if _internal_max_reps is None:
        _internal_max_reps = n_reps
    elif _internal_max_reps > n_reps:
        raise CheckError(
            f"Specified _internal_max_reps to be {_internal_max_reps} "
            f"however each fit only has {n_reps} replicas"
        )

    if _internal_min_reps is None:
        _internal_min_reps = 10

    if _internal_max_reps < _internal_min_reps:
        raise CheckError(
            f"maximum replicas per fit, {_internal_max_reps}, is less than minimum replicas "
            f", {_internal_min_reps}. If you have set _internal_max_reps and"
            "_internal_min_reps then ensure that they take sensible values."
        )
    return {
        "_internal_max_reps": _internal_max_reps,
        "_internal_min_reps": _internal_min_reps,
    }


@make_argcheck
def check_fits_different_filterseed(fits):
    """Input fits should have the different filter seed if they are being
    used for multiple closure test studies, because in high-level hands-waving
    terms the different level 1 shifts represents different
    'runs of the universe'!

    """
    seed_list = [fit.as_input()["closuretest"]["filterseed"] for fit in fits]
    if len(seed_list) > len(set(seed_list)):
        raise CheckError(
            f"Multiclosure actions require that fits had different level 1 noise"
        )
