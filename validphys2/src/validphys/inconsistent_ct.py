import dataclasses
from validphys.coredata import CommonData
import pandas as pd



@dataclasses.dataclass(eq = False)
class InconsistentCommonData(CommonData):
    """
    Class that inherits all of the methods
    of coredata.CommonData class.
    
    This class is meant to have all the 
    methods needed in order to introduce 
    an inconsistency within a Closure Test.
    """

    setname: str
    ndata: int
    commondataproc: str
    nkin: int
    nsys: int
    commondata_table: pd.DataFrame = dataclasses.field(repr=False)
    systype_table: pd.DataFrame = dataclasses.field(repr=False)
    systematics_table: pd.DataFrame = dataclasses.field(init=None, repr=False)


    def with_MULT_sys(self,mult_sys):
        """
        returns an InconsistentCommonData instance 
        with MULT systematics replaced by mult_sys

        Parameters
        ----------
        mult_sys : pd.DataFrame()
                 all MULT columns of 
                 InconsistentCommonData.commondata_table
        """
        table = self.commondata_table.copy()
        table["MULT"] = mult_sys
        return dataclasses.replace(self, commondata_table = table)
    
    def with_ADD_sys(self,add_sys):
        """
        returns an InconsistentCommonData instance
        with ADD systematics replaced by add_sys

        Parameters
        ----------
        add_sys : pd.DataFrame()
                 all ADD columns of
                 InconsistentCommonData.commondata_table
        """
        table = self.commondata_table.copy()
        table["ADD"] = add_sys
        return dataclasses.replace(self, commondata_table = table)

    def rescale_sys(self,type_err,CORR,UNCORR,sys_rescaling_factor):
        """
        rescale the sys (MULT or ADD) by constant factor, sys_rescaling_factor,
        a distinction is done between CORR and UNCORR systematics

        Parameters
        ----------

        type_err : str 
                e.g. 'MULT' or 'ADD'

        CORR : bool

        UNCORR : bool

        sys_rescaling_factor : float, int

        Returns
        -------
        pd.DataFrame corresponding to the rescaled MULT systematics
        """
        
        err_table = self.systematics_table.loc[:,[type_err]].copy()
        # get indices of CORR / UNCORR sys
        systype_corr = self.systype_table[(self.systype_table["type"] == type_err) 
                            & (~self.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]
        
        systype_uncorr = self.systype_table[(self.systype_table["type"] == type_err) 
                            & (self.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]

        # rescale systematics
        if CORR:
            err_table.iloc[:,systype_corr.index - 1] *= sys_rescaling_factor
        if UNCORR:
            err_table.iloc[:,systype_uncorr.index - 1] *= sys_rescaling_factor

        return err_table
    
    def process_commondata(
                    self,ADD,MULT,
                    CORR,UNCORR,
                    inconsistent_datasets,
                    sys_rescaling_factor
                    ):
        """
        returns a commondata instance
        with modified systematics. 
        Note that if commondata.setname
        is not within the inconsistent_datasets or if both ADD and
        MULT are False, then the commondata object
        will not be modified.

        Parameters
        ----------

        ADD : bool

        MULT : bool

        CORR : bool

        UNCORR : bool

        inconsistent_datasets : list
                            list of the datasets for which an inconsistency should be introduced
        
        sys_rescaling_factor : float, int

        Returns
        -------
        validphys.inconsistent_ct.InconsistentCommonData
        """
        new_commondata = self

        if not self.setname in inconsistent_datasets:
            return self
        
        # When modifying new_commondata I think one should call the method with_MULT/ADD on the new_commondata
        # itself. This for the fact that if one has both MULT and ADD set to TRUE the modified commondata would
        # takes into account only the ADD keyword
        if MULT:
            new_commondata = new_commondata.with_MULT_sys(self.rescale_sys("MULT",CORR,UNCORR,sys_rescaling_factor))
            
        if ADD:
            new_commondata = new_commondata.with_ADD_sys(self.rescale_sys("ADD",CORR,UNCORR,sys_rescaling_factor))

        return new_commondata

