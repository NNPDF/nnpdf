from validphys.commondataparser import parse_commondata_new
from validphys.api import API
from validphys.covmats import covmat_from_systematics 

import numpy as np

datasetname = "CMS_2JET_7TEV_M12Y"

cd_new = parse_commondata_new(datasetname, variants=[])

inp = {'dataset_input': {"dataset": "CMS_2JET_7TEV"} , 'theoryid': 200, 'use_cuts': 'internal'}

ds = API.dataset(**inp)
cd_old = ds.load_commondata()

covmat_new = covmat_from_systematics(cd_new, ds, use_weights_in_covmat=False)

covmat_old = covmat_from_systematics(cd_old, ds, use_weights_in_covmat=False)

print("New Covariance Matrix is the same as Old Covariance Matrix:")
print(np.allclose(covmat_new, covmat_old))
print()
print(covmat_old / covmat_new)
print(cd_old)

print(cd_new)

