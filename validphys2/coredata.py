"""
Data containers backed by Python managed memory (Numpy arrays and Pandas
dataframes).  This module is intended to substitute large parts of the C++
wrappers.
"""
import dataclasses
import numpy as np
import pandas as pd

@dataclasses.dataclass(eq=False)
class CommonData:
    """
    Data contained in Commondata files
    Parameters
    ----------
    setname: str
        Name of the dataset
    ndata: int
        Number of data points
    data: array of floats with length ndata
        Data values
    commondataproc: str
        Process type, one of 21 options. 
    nkin: int
        Number of kinematics specified
    kinematics: list of str with length nkin
        Kinematic variables kin1, kin2, kin3 ...
    nsys: int
        Number of systematics
    sysid: list of str with length nsys
        ID for systematic
    stat: array of floats with length ndata
        Statistical uncertainties on each data point
        (separate ADD and MULT here?)
    sys: array of floats with dimensions ndat x nsys
        Systematic uncertainties on each data point
        (separate ADD and MULT here?)
    

    
    """
    setname: str
    ndata: int
    data: np.array 
    commondataproc: str
    nkin: int 
    kinematics: list(str)
    nsys: int
    sysid: list(str)
    stat: np.array
    sys: np.array 
