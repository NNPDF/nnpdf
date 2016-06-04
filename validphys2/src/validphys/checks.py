# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 19:35:40 2016

@author: Zahari Kassabov
"""
import tempfile

import lhapdf

from reportengine.checks import make_check, CheckError, require_one, check_not_empty

@make_check
def check_pdf_is_montecarlo(ns, **kwargs):
    pdf = ns['pdf']
    etype = pdf.ErrorType
    if etype != 'replicas':
        raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
                          % (pdf, etype))

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
