from validphys.commondataparser import parse_commondata_new
from validphys.api import API
from validphys.covmats import covmat_from_systematics 

import numpy as np

datasetname = "ATLAS_1JET_8TEV_R06_PTY"

cd_new = parse_commondata_new(datasetname, variants=[])

inp = {'dataset_input': {"dataset": "ATLAS_1JET_8TEV_R06"} , 'theoryid': 200, 'use_cuts': 'internal'}

ds = API.dataset(**inp)
cd_old = ds.load_commondata()

covmat_new = covmat_from_systematics(cd_new, ds, use_weights_in_covmat=False)

covmat_old = covmat_from_systematics(cd_old, ds, use_weights_in_covmat=False)

print(np.allclose(covmat_new, covmat_old))
print(cd_old)

print(cd_new)

