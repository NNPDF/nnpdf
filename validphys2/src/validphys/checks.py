# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 19:35:40 2016

@author: Zahari Kassabov
"""
from reportengine.checks import make_check, CheckError, require_one, check_not_empty

@make_check
def check_pdf_is_montecarlo(ns, **kwarks):
    pdf = ns['pdf']
    etype = pdf.ErrorType
    if etype != 'replicas':
        raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
                          % (pdf, etype))
