# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 19:35:40 2016

@author: Zahari Kassabov
"""
import tempfile

from matplotlib import scale as mscale

import lhapdf

from reportengine.checks import (make_check, CheckError, require_one,
                                 check_not_empty, make_argcheck)

from validphys import lhaindex

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

#TODO: Make postfit know about which replicas it has selected
@make_check
def check_has_fitted_replicas(ns, **kwargs):
    name, path = ns['fit']
    postfit_path = path/'nnfit'/'postfit.log'
    if not postfit_path.exists():
        raise CheckError("Fit {name} does not appear to be completed. "
        "Expected to find file {postfit_path}".format(**locals()))

    if not lhaindex.isinstalled(name):
        raise CheckError("The PDF corresponding to the fit, '%s', "
        "needs to be "
        "installed in LHAPDF (i.e. copied to %s)."%
        (name, lhaindex.get_lha_datapath()))


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

