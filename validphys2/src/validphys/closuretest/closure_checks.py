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
