

def process_commondata(commondata,ADD,MULT,CORR,UNCORR,inconsistent_datasets,sys_rescaling_factor):
    """
    Given a commondata instance return a commondata instance
    with modified systematics. Note that if commondata.setname
    is not within the inconsistent_datasets or if both ADD and
    MULT are False, then the input commondata is outputted

    Parameters
    ----------
    
    commondata : validphys.coredata.CommonData

    ADD : bool

    MULT : bool

    CORR : bool

    UNCORR : bool

    inconsistent_datasets : list
                        list of the datasets for which an inconsistency should be introduced
    
    sys_rescaling_factor : float, int

    Returns
    -------
    validphys.coredata.CommonData
    """

    if not commondata.setname in inconsistent_datasets:
        return commondata
    
    if MULT and ADD:
        cd_tmp = commondata.with_MULT_sys(commondata.multiplicative_errors_rescale(CORR,UNCORR,sys_rescaling_factor))
        cd = cd_tmp.with_ADD_sys(cd_tmp.additive_errors_rescale(CORR,UNCORR,sys_rescaling_factor))
        return cd
    
    if MULT:
        cd = commondata.with_MULT_sys(commondata.multiplicative_errors_rescale(CORR,UNCORR,sys_rescaling_factor))
        return cd
   
    if ADD:
        cd = commondata.with_ADD_sys(commondata.additive_errors_rescale(CORR,UNCORR,sys_rescaling_factor))
        return cd
    
    return commondata