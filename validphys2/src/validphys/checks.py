# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 19:35:40 2016

@author: Zahari Kassabov
"""
import tempfile
from collections import Counter

from matplotlib import scale as mscale

import lhapdf

from reportengine.checks import (make_check, CheckError, require_one,
                                 check_not_empty, make_argcheck, check_positive, check)

from validphys import lhaindex
from validphys.core import CutsPolicy

import logging
log = logging.getLogger(__name__)

@make_check
def check_pdf_is_montecarlo(ns, **kwargs):
    pdf = ns['pdf']
    etype = pdf.ErrorType
    if etype != 'replicas':
        raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
                          % (pdf, etype))

@make_check
def check_know_errors(ns, **kwargs):
    pdf = ns['pdf']
    try:
        pdf.nnpdf_error
    except NotImplementedError as e:
        raise CheckError(e) from e


@make_check
def check_can_save_grid(ns, **kwags):
    if not ns['installgrid']:
        return

    write_path = lhapdf.paths()[-1]
    try:
        tempfile.TemporaryFile(dir=write_path)
    except OSError as e:
        raise CheckError("Cannot write to the LHAPDF path %s.\n"
                         "This is required because the 'installgrid' "
                         "parameter is set to True:\n%s" %
                        (write_path, e))

@make_argcheck
def check_xlimits(xmax, xmin):
    if not (0 <= xmin < xmax <= 1):
        raise CheckError(f'xmin ({xmin}) and xmax ({xmax}) must satisfy \n'
                         '0 <= xmin < xmax <= 1')

@make_check
def check_has_fitted_replicas(ns, **kwargs):
    name, path = ns['fit']
    postfit_path = path/'postfit'/'postfit.log'
    old_postfit_path = path/'nnfit'/'postfit.log'
    if not postfit_path.exists():
        if not old_postfit_path.exists():
            raise CheckError(
                f"Fit {name} does not appear to be completed. "
                f"Expected to find file {postfit_path}")
        else:
            log.info(f"Cannot find postfit log at: {postfit_path}. "
                     f"Falling back to old location: {old_postfit_path}")

    if not lhaindex.isinstalled(name):
        raise CheckError(
            f"The PDF corresponding to the fit, '{name}'"
            "needs to be installed in LHAPDF (i.e. copied to "
            f"{lhaindex.get_lha_datapath()})."
        )


def check_scale(scalename, allow_none=False):
    """Check that we have a valid matplotlib scale. With allow_none=True,
    also None is valid."""
    @make_check
    def check(ns, *args, **kwargs):
        val = ns[scalename]
        if val is None and allow_none:
            return
        valid_scales = mscale.get_scale_names()
        if not val in valid_scales:
            e = CheckError("Invalid plotting scale: %s" % scalename,
                             bad_item=val,
                             alternatives=valid_scales,
                             display_alternatives='all')
            e.alternatives_header = "No such scale '%s'. The allowed values are:"
            raise e
    return check

@make_argcheck
def check_cuts_considered(use_cuts):
    if use_cuts not in (CutsPolicy.FROMFIT, CutsPolicy.INTERNAL):
        raise CheckError(f"Cuts must be computed for this action, but they are set to {use_cuts.value}")


@make_argcheck
def check_have_two_pdfs(pdfs):
    check(len(pdfs) == 2,'Expecting exactly two pdfs.')


#The indexing to one instead of zero is so that we can be consistent with
#how plot_fancy works, so normalize_to: 1 would normalize to the first pdf
#for both.
@make_argcheck
def check_pdf_normalize_to(pdfs, normalize_to):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be, a pdf id or an index of the "
           "pdf (starting from one)")

    if normalize_to is None:
        return

    names = [pdf.name for pdf in pdfs]
    if isinstance(normalize_to, int):
        normalize_to -= 1
        if not normalize_to < len(names) or normalize_to<0:
            raise CheckError(msg)
        return {'normalize_to': normalize_to}

    if isinstance(normalize_to, str):
        try:
            normalize_to = names.index(normalize_to)
        except ValueError:
            raise CheckError(msg, normalize_to, alternatives=names)
        return {'normalize_to': normalize_to}


    raise RuntimeError("Should not be here")


def _check_list_different(l, name):
    strs = [str(item) for item in l]
    if not len(set(strs))==len(l):
        counter = Counter(strs)
        duplicates = [k for k, v in counter.items() if v > 1]
        raise CheckError(f"{name} must be all different "
                         f"but there are duplicates: {duplicates}")

@make_argcheck
def check_fits_different(fits):
    """Need this check because oterwise the pandas object gets confused"""
    return _check_list_different(fits, 'fits')

@make_argcheck
def check_dataspecs_fits_different(dataspecs_fit):
    """Need this check because oterwise the pandas object gets confused"""
    return _check_list_different(dataspecs_fit, 'fits')

@make_argcheck
def check_speclabels_different(dataspecs_speclabel):
    """This is needed for grouping dataframes (and because
    generally indecated a bug)"""
    return _check_list_different(dataspecs_speclabel, 'dataspecs_speclabel')

@make_argcheck
def check_two_dataspecs(dataspecs):
    l = len(dataspecs)
    check(l == 2, f"Expecting exactly 2 dataspecs, not {l}")
