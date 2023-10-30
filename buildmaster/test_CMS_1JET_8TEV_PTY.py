from validphys.commondataparser import parse_new_metadata, parse_commondata_new
from validphys.api import API
from validphys.covmats import covmat_from_systematics 

import numpy as np
import pathlib

datasetname = "CMS_1JET_8TEV_PTY"
metadata_file = pathlib.Path("/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/CMS_1JET_8TEV/metadata.yaml")
observable_name = "PTY"
metadata = parse_new_metadata(metadata_file, observable_name)

# cd_new = parse_commondata_new(datasetname, variants=[])
cd_new = parse_commondata_new(metadata)
#inp = {'dataset_input': {"dataset": "CMS_1JET_8TEV"} , 'theoryid': 200, 'use_cuts': 'internal'}

inps = [{'dataset': "CMS_1JET_8TEV"}]
inp = dict(dataset_inputs=inps, theoryid=200, use_cuts="internal")
covmat_old = API.dataset_inputs_covmat_from_systematics(**inp)

inp = {'dataset_input': {"dataset": "CMS_1JET_8TEV"} , 'theoryid': 200, 'use_cuts': 'internal'}

ds = API.dataset(**inp)
#cd_old = ds.load_commondata()

covmat_new = covmat_from_systematics(cd_new, ds, use_weights_in_covmat=False)

#covmat_old = covmat_from_systematics(cd_old, ds, use_weights_in_covmat=False)
ones = covmat_new / covmat_old

print("New Covariance Matrix is the same as Old Covariance Matrix:")
print(np.allclose(ones, np.ones(covmat_old.shape)))
print()
print(ones)
print()
print(np.diag(ones))
print()
print(np.min(ones), np.max(ones))

#print(cd_old)
#
#print(cd_new)
#
