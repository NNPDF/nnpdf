"""
closuretest/checks.py

Module containing checks specific to the closure tests.

"""
from collections import defaultdict
import logging

from reportengine.checks import CheckError, make_argcheck

log = logging.getLogger(__name__)


@make_argcheck
def check_use_fitcommondata(use_fitcommondata):
    """Base check that `use_fitcommondata` is being used, check should be used
    with all actions which require comparison to fitcommondata
    """
    if not use_fitcommondata:
        raise CheckError("use_fitcommondata must be set to True for closure test estimators")


@make_argcheck
def check_fit_isclosure(fit):
    """Check the input fit is a closure test"""
    fitinfo = fit.as_input()
    if not "closuretest" in fitinfo:
        raise CheckError(
            f"There is no `closuretest` namespace in {fit}'s runcard. "
            f"{fit} is therefore not suitable for closure-test studies."
        )
    if not "fakedata" in fitinfo["closuretest"]:
        raise CheckError(
            f"The `fakedata` key does not exist in the `closuretest` namespace of {fit}'s runcard. "
            f"{fit} is therefore not suitable for closure-test studies."
        )
    if not fitinfo["closuretest"]["fakedata"]:
        raise CheckError(
            f"The `fakedata` key is not set to `true` in the `closuretest` namespace of {fit}'s runcard. "
            f"{fit} is therefore not suitable for closure-test studies."
        )


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
        check_fit_isclosure.__wrapped__(fit)


@make_argcheck
def check_t0pdfset_matches_law(t0pdfset, fit):
    t0_from_fit = fit.as_input()["closuretest"]["fakepdf"]
    if not str(t0pdfset) == t0_from_fit:
        raise CheckError(f"Underlying pdf: {t0_from_fit}, does not match t0pdfset: {t0pdfset}")


@make_argcheck
def check_t0pdfset_matches_multiclosure_law(multiclosure_underlyinglaw, t0set):
    """Checks that, if a multiclosure_underlyinglaw is present, it matches the t0set
    Checks t0set instead of t0pdfset since different mechanisms can fill t0set
    """
    if str(t0set) != str(multiclosure_underlyinglaw):
        log.warning(
            f"The underlying pdf {multiclosure_underlyinglaw} does not match t0pdfset: {t0set}"
        )


@make_argcheck
def check_at_least_10_fits(fits):
    if len(fits) < 10:
        raise CheckError(
            "Multiclosure actions testing finite sampling effects require at least 10 fits"
        )


@make_argcheck
def check_multifit_replicas(fits_pdf, _internal_max_reps, _internal_min_reps):
    """Checks that all the fit pdfs have the same number of replicas N_rep. Then
    check that N_rep is greater than the smallest number of replicas used in
    actions which subsample the replicas of each fit.

    This check also has the secondary
    effect of filling in the namespace key _internal_max_reps
    which can be used to override the number of replicas used at the level of
    the runcard, but by default get filled in as the number of replicas in each
    fit.

    """
    n_reps = {pdf.get_members() - 1 for pdf in fits_pdf}
    if len(n_reps) != 1:
        raise CheckError("all fits for multiclosure actions should have same number of replicas")
    n_reps = n_reps.pop()
    if _internal_max_reps is None:
        _internal_max_reps = n_reps
    elif _internal_max_reps > n_reps:
        raise CheckError(
            f"Specified _internal_max_reps to be {_internal_max_reps} "
            f"however each fit only has {n_reps} replicas"
        )

    if _internal_max_reps < _internal_min_reps:
        raise CheckError(
            f"maximum replicas per fit, {_internal_max_reps}, is less than minimum replicas "
            f", {_internal_min_reps}. If you have set _internal_max_reps and"
            "_internal_min_reps then ensure that they take sensible values."
        )
    return {"_internal_max_reps": _internal_max_reps, "_internal_min_reps": _internal_min_reps}


@make_argcheck
def check_fits_different_filterseed(fits):
    """Input fits should have different filter seeds if they are being
    used for multiple closure test studies, because in high-level hand-waving
    terms the different level 1 shifts represents different
    'runs of the universe'!

    """
    seed_fits_dict = defaultdict(list)

    for fit in fits:
        pdf_name = fit.as_input()["pdf"]["id"]
        seed = fit.as_input()["closuretest"]["filterseed"]
        seed_fits_dict[seed].append(pdf_name)

    bad_fits = [fits for _, fits in seed_fits_dict.items() if len(fits) > 1]

    if bad_fits:
        raise CheckError(
            "Multiclosure actions require that fits have different level 1 "
            "noise and therefore different filter seeds. The following groups "
            f"of fits have the same seed: {bad_fits}."
        )


@make_argcheck
def check_fits_have_same_basis(fits_basis):
    """Check the basis is the same for all fits"""
    bases_set = set(fits_basis)
    if len(bases_set) != 1:
        raise CheckError(
            "All fits must have the same basis, however the following bases "
            f"were found: {bases_set}"
        )
