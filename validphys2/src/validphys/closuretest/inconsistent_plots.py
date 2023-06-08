"""
closuretest/inconsistent_plots.py

Useful plots for analysis of inconsistent closure tests
"""
from reportengine.figure import figure
from reportengine.table import table
from validphys import plotutils
from validphys.inconsistent_ct import InconsistentCommonData
from validphys.covmats import dataset_inputs_covmat_from_systematics
import numpy as np
import pandas as pd

def covmat_trace(dataset_input_list,commondata_wc):
    """Return trace of experimental matrix
    """
    normal_covmat = dataset_inputs_covmat_from_systematics(
                commondata_wc,
                dataset_input_list,
                use_weights_in_covmat=False,
                norm_threshold=None,
                _list_of_central_values=None,
                _only_additive=False,
            )
    return np.trace(normal_covmat)

def mod_covmat_trace(dataset_input_list,commondata_wc, inconsistent_datasets, ADD, MULT, CORR, UNCORR, SPECIAL, sys_rescaling_factor):
    """ Calculate trace of inconsistent covmat rescaled by
    sys_rescaling_factor affecting ADD/MULT & CORR/UNCORR/SPECIAL.
    """
    commondata_wc_temp = [cd.process_commondata(ADD,MULT,CORR,UNCORR,SPECIAL,
                                            inconsistent_datasets,sys_rescaling_factor)
                            for cd in commondata_wc]
    modified_covmat = dataset_inputs_covmat_from_systematics(
            commondata_wc_temp,
            dataset_input_list,
            use_weights_in_covmat=False,
            norm_threshold=None,
            _list_of_central_values=None,
            _only_additive=False,
        )
    #Calculate trace of modified trace (either for type 1 or 2)
    modified_trace = np.trace(modified_covmat)
    return modified_trace

@figure
def plot_trace_impact(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL):
    """
    Plot trace ratio for different sys_rescaling_factors. Specify what kind 
    of error has been modified and for which datasets.
    The marked points are the one for which the trace ratio corresponds between type1/type2 inconsistent fit.
    """

    # Load here all the data, does not make sense to load them each time the funciton is called
    dataset_input_list = list(data.dsinputs)
    commondata_wc = data.load_commondata_instance()
    commondata_wc_ic = [
                InconsistentCommonData(setname=cd.setname, ndata=cd.ndata, 
                    commondataproc=cd.commondataproc, 
                    nkin=cd.nkin, nsys=cd.nsys, 
                    commondata_table = cd.commondata_table, 
                    systype_table = cd.systype_table) 
                for cd in commondata_wc
                ]
    normal_trace = covmat_trace(dataset_input_list,commondata_wc)
    lam_factors = np.arange(0,3,0.02)
    ratios = []
    fig, ax = plotutils.subplots()
    points = []
    i = 0
    for lam in lam_factors:
        mod_trace = mod_covmat_trace(dataset_input_list,commondata_wc_ic, inconsistent_datasets, 
                                     ADD,MULT,CORR,UNCORR,SPECIAL, 
                                     lam)
        if lam < 1: ratios.append(mod_trace/normal_trace*100)
        if lam >= 1: ratios.append(normal_trace/mod_trace*100)
        if i%10 == 0 and lam < 1: 
            ax.plot(lam_factors[i],ratios[-1],marker = ".",
                    markersize = 10, 
                    label = "lambda type 2: " + str(round(lam_factors[i],3)) + "; ratio: " + str(round(ratios[-1],3)))
            points.append(ratios[-1])
        i += 1
    for point in points:
        a,b = find_intersections(np.asarray(lam_factors), np.asarray(ratios), point)
        ax.plot(a,b,marker = ".",markersize = 10,label = "lambda type1: " + str(round(a[0],3))+"; ratio: " + str(round(b[0],3)))
    type_a_m = ""
    type_c_u_s = ""
    if ADD: type_a_m = " ADD "
    if MULT: type_a_m = type_a_m + " MULT "
    if CORR: type_c_u_s = " CORR "
    if UNCORR: type_c_u_s = type_c_u_s + " UNCORR "
    if SPECIAL: type_c_u_s = type_c_u_s  + " SPECIAL"
    ax.plot(lam_factors, ratios, label = "percentage ratios")
    title = "Impact of inconsistency of type " + str(type_a_m)  + " and " + str(type_c_u_s) + " in \n" + str(inconsistent_datasets) + " wrt all ds. \n"
    ax.legend()
    ax.set_title(title)
    ax.set_xlabel("rescaling factor")
    ax.set_ylabel("percentage ratio")
    return fig

def find_intersections(x, y, C):
    # Contains numpy indexing tricks that can be hard to reproduce
    # in the case where both functions are non-constants
    ii, = np.nonzero((y[1:]-C)*(y[:-1]-C) < 0.)  # intersection indices
    x_intersections = x[ii] + (C - y[ii])/(y[1+ii] - y[ii])*(x[1+ii] - x[ii])
    y_intersections = C * np.ones(len(ii))
    return x_intersections, y_intersections