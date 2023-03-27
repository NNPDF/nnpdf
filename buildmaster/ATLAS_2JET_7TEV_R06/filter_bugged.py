from validphys.loader import Loader
import numpy as np
import yaml


INTRA_DATASET_SYS_NAME = ("UNCORR", "CORR", "THEORYUNCORR", "THEORYCORR")
setname = "ATLAS_2JET_7TEV_R06"

def decompose_covmat(covmat):
     """Given a covmat it return an array sys with shape (ndat,ndat)
     giving ndat correlated systematics for each of the ndat point.
     The original covmat is obtained by doing sys@sys.T"""

     lamb, mat = np.linalg.eig(covmat)
     sys = np.multiply(np.sqrt(lamb), mat)
     return sys


def filter_ATLAS_2JET_7TEV_R06_uncertainties_bugged():
    """
    - read covmat from old CommonData format
    - compute artificial uncertainties and save them to
      .yaml file
    """
    
    l = Loader() 
    cd = l.check_commondata(setname=setname).load_commondata_instance()

    sys_errors = cd.systematic_errors()
    
    is_intra_dataset_error = sys_errors.columns.isin(INTRA_DATASET_SYS_NAME)

    sys = sys_errors.loc[:,is_intra_dataset_error]

    # construct sys Covmat and decompose so as to have the same
    # amount of artificial sys as in the new CD implementation
    C_sys = np.einsum('ij,kj->ik',sys,sys)
    A_sys = decompose_covmat(C_sys)
    
    # special sys, ATLASLUMI11
    atlas_lumi = sys_errors.loc[:, ~is_intra_dataset_error].to_numpy()

    # error definition
    error_definition = {f"art_sys_{i}" : {"description":  f"artificial systematic {i}",
                                        "treatment": "ADD", "type": "CORR"}
                        for i in range(1,A_sys.shape[0]+1)}

    error_definition["luminosity_uncertainty"] = {"description": "luminosity uncertainty",
                                                    "treatment": "ADD", "type": "ATLASLUMI11"}

    # store error in dict
    error = []
    for n in range(A_sys.shape[0]):
        error_value={}
        for m in range(A_sys.shape[1]):
            error_value[f"art_sys_{m+1}"] = float(A_sys[n,m])
        
        error_value["luminosity_uncertainty"] = float(atlas_lumi[n])
        error.append(error_value)

    uncertainties_yaml = {"definition": error_definition, "bins": error}
    
    # save uncertainties to .yaml file
    with open("uncertainties_bugged.yaml",'w') as file:
        yaml.dump(uncertainties_yaml,file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_2JET_7TEV_R06_uncertainties_bugged()