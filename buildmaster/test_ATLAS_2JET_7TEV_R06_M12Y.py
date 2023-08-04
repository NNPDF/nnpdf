from validphys.commondataparser import parse_commondata_new
from validphys.api import API
from validphys import covmats 

import numpy as np

datasetname = "ATLAS_2JET_7TEV_R06_M12Y"
inp = {'dataset_input': {"dataset": "ATLAS_2JET_7TEV_R06"} , 'theoryid': 200, 'use_cuts': 'internal'}
ds = API.dataset(**inp)

cd_new = parse_commondata_new(datasetname, variants=['bugged'])

cd_old = ds.load_commondata()

covmat_new = covmats.covmat_from_systematics(cd_new, ds, use_weights_in_covmat=False)
sqrt_covmat_new = covmats.sqrt_covmat(covmat_new)

covmat_old = covmats.covmat_from_systematics(cd_old, ds, use_weights_in_covmat=False)
sqrt_covmat_old = covmats.sqrt_covmat(covmat_old)


print("New Covariance Matrix is the same as Old Covariance Matrix:")
print(np.allclose(covmat_new, covmat_old))
print()
print(covmat_old / covmat_new)
print(cd_old)

print(cd_new)

