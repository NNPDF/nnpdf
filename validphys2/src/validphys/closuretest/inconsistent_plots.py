"""
closuretest/inconsistent_plots.py

Useful plots for analysis of inconsistent closure tests
"""
from reportengine.figure import figure
from validphys import plotutils
from validphys.inconsistent_ct import InconsistentCommonData
from validphys.covmats import dataset_inputs_covmat_from_systematics
import numpy as np

def covmat_trace_ratio(data, inconsistent_datasets, ADD, MULT, CORR, UNCORR, SPECIAL, sys_rescaling_factor, type1 = True, type2 = False):
    """ Calculate trace difference between the original exp covmat and the inconsistent one rescaled by
    sys_rescaling_factor affecting ADD/MULT (decided by type_a_m) unc CORR/UNCORR/SPECIAL (decided by type_c_u_s)

    Separate two possible cases: type1 and type2 trace diff.

    This is necessary in order to have a good comparison useful for the inconsistent closure tests of type1 and
    type2. 

    A type1 reference fit is defined by: 
    L1_data = L0 + N(0,C+)
    L2_data = L1 + N(0,C+)

    A type2 reference fit is defined by:

    L1_data = L0 + N(0,C_exp)
    L2_data = L1 + N(0,C_exp)

    where C+ is (roughly) defined as C_exp*sys_rescaling_factor (usually larger than 1)

    in order to see the impact of the inconsistency then:

    Type1 inconsistent closure test:
    L1_data = L0 + N(0,C+)
    L2_data = L1 + N(0,C_exp)

    Type1 trace percentage diff: 

    tr(C_exp)/(tr(C+))*100

    Type2 inconsistent closure test:
    L1_data = L0 + N(0,C_exp)
    L2_data = L1 + N(0,C-)

    where C- is (roughly) defined as C_exp*sys_rescaling_factor (usually smaller than 1)


    Type2 percentage diff:
    tr(C-)/(tr(C_exp))*100


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
    modified_trace = np.trace(modified_covmat)
    if type1:   
        percentage_ratio = (trace)/(modified_trace)*100
    if type2:
        percentage_ratio = modified_trace/trace*100
    return percentage_ratio

@figure
def plot_trace_impact(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL):
    rs_factors_2 = np.arange(0.1,1.1,0.1)
    rs_factors_1 = np.arange(1,5,0.1)
    impacts_1 = []
    impacts_2 = []
    for fac in rs_factors_2:
        impacts_2.append(covmat_trace_ratio(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL, fac,type1 = False,type2 = True))
    for fac in rs_factors_1:
        impacts_1.append(covmat_trace_ratio(data, inconsistent_datasets, ADD,MULT,CORR,UNCORR,SPECIAL, fac,type1 = True,type2 = False))
    fig, ax = plotutils.subplots()

    ax.plot(rs_factors_1,impacts_1, label = "type1 Inconsistent Closure")
    ax.plot(rs_factors_2,impacts_2, label = "type2 Inconsistent Closure")
    type_a_m = ""
    type_c_u_s = ""
    if ADD: type_a_m = " ADD "
    if MULT: type_a_m = type_a_m + " MULT "
    if CORR: type_c_u_s = " CORR "
    if UNCORR: type_c_u_s = type_c_u_s + " UNCORR "
    if SPECIAL: type_c_u_s = type_c_u_s  + " SPECIAL"
    title = "Impact of inconsistency of type " + str(type_a_m)  + " and " + str(type_c_u_s) + " in \n" + str(inconsistent_datasets) + " wrt all ds. \n"
    ax.legend()
    ax.set_title(title)
    ax.set_xlabel("rescaling factor")
    ax.set_ylabel("percentage ratio")
    return fig