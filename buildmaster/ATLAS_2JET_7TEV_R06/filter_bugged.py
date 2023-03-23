from validphys.api import API
import numpy as np
import yaml

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
    inp = dict(
     dataset_input={'dataset': 'ATLAS_2JET_7TEV_R06'},
     theoryid=200,
     use_cuts="internal"
     )

    covmat = API.covmat_from_systematics(**inp)

    # generate artificial systematics 
    A_art_sys = decompose_covmat(covmat=covmat)
    
    # store error in dict
    error = []
    for n in range(A_art_sys.shape[0]):
        error_value={}
        for m in range(A_art_sys.shape[1]):
            error_value[f"art_sys_{m+1}"] = float(A_art_sys[n,m])
        error.append(error_value)

    # error definition
    error_definition = {f"art_sys_{i}" : {"description":  f"bugged artificial systematic {i}",
                                        "treatment": "ADD", "type": "CORR"}
                        for i in range(1,A_art_sys.shape[0]+1)}
    
    uncertainties_yaml = {"definition": error_definition, "bins": error}
    
    
    with open("uncertainties_bugged.yaml",'w') as file:
        yaml.dump(uncertainties_yaml,file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_2JET_7TEV_R06_uncertainties_bugged()