"""
closuretest/inconsistent_plots.py

Useful plots for analysis of inconsistent closure tests
"""
from reportengine.figure import figure
from validphys import plotutils
from validphys.inconsistent_ct import InconsistentCommonData
from validphys.covmats import dataset_inputs_covmat_from_systematics
import numpy as np

#TODO: maybe change location
def covmat_trace_diff(data, inconsistent_datasets, ADD, MULT, CORR, UNCORR, SPECIAL, sys_rescaling_factor):
    """ Calculate trace difference between the original exp covmat and the inconsistent one rescaled by
    sys_rescaling_factor affecting ADD/MULT (decided by type_a_m) unc CORR/UNCORR/SPECIAL (decided by type_c_u_s)
    """

    dataset_input_list = list(data.dsinputs)
    commondata_wc = data.load_commondata_instance()
    commondata_wc = [
                    InconsistentCommonData(setname=cd.setname, ndata=cd.ndata, 
                        commondataproc=cd.commondataproc, 
                        nkin=cd.nkin, nsys=cd.nsys, 
                        commondata_table = cd.commondata_table, 
                        systype_table = cd.systype_table) 
                    for cd in commondata_wc
                    ]
    consistent_covmat = dataset_inputs_covmat_from_systematics(
        commondata_wc,
        dataset_input_list,
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
    )
    
    trace = np.trace(consistent_covmat)
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
    inconsistent_trace = np.trace(modified_covmat)
    percentage_diff = abs(inconsistent_trace-trace)/(trace)*100
    return percentage_diff

@figure
def plot_trace_impact(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL):
    rs_factors = np.arange(0.1,5,0.1)
    impacts = []
    for fac in rs_factors:
        impacts.append(covmat_trace_diff(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL, fac))
    fig, ax = plotutils.subplots()

    ax.plot(rs_factors,impacts)
    ds_names = ""
    for ds in inconsistent_datasets:
        ds_names += str(ds) + " "
    type_a_m = ""
    type_c_u_s = ""
    if ADD: type_a_m = " ADD "
    if MULT: type_a_m = type_a_m + " MULT "
    if CORR: type_c_u_s = " CORR "
    if UNCORR: type_c_u_s = type_c_u_s + " UNCORR "
    if SPECIAL: type_c_u_s = type_c_u_s  + " SPECIA L"
    title = "Impact of inconsistency of type " + str(type_a_m)  + " and " + str(type_c_u_s) + " in " + ds_names + " wrt all ds"
    ax.set_title(title)
    ax.set_xlabel("rescaling factor")
    ax.set_ylabel("percentage diff")
    return fig
