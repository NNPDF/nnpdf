from validphys.commondataparser import parse_commondata_new
from validphys.api import API
from validphys import covmats
from validphys.calcutils import calc_chi2

import numpy as np

datasetname = "ATLAS_1JET_8TEV_R06_PTY"

inp = {
    'dataset_input': {"dataset": "ATLAS_1JET_8TEV_R06_DEC"},
    'theoryid': 200,
    'use_cuts': 'internal',
    'pdf': "NNPDF40_nnlo_as_01180",
}

ds = API.dataset(**inp)
pdf = API.pdf(**inp)

cd_new = parse_commondata_new(datasetname, variants=['decorrelated'])

cd_old = ds.load_commondata()

covmat_new = covmats.covmat_from_systematics(cd_new, ds, use_weights_in_covmat=False)
sqrt_covmat_new = covmats.sqrt_covmat(covmat_new)

covmat_old = covmats.covmat_from_systematics(cd_old, ds, use_weights_in_covmat=False)
sqrt_covmat_old = covmats.sqrt_covmat(covmat_old)

diff = cd_old.central_values - covmats.dataset_t0_predictions(ds, pdf)

print(f"Old Covmat == New Covmat : {np.allclose(covmat_new, covmat_old)}")
print()
print(
    f"Old Central Values == New Central Values: {np.allclose(cd_old.central_values, cd_new.central_values)}"
)
print()
print(
    f"Old Chi2 {calc_chi2(sqrt_covmat_old, diff) / cd_old.ndata:.3f}, New Chi2 {calc_chi2(sqrt_covmat_new, diff)/cd_new.ndata:.3f}"
)
print()
import matplotlib.pyplot as plt
import seaborn as sns

fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# Create a shared colorbar axis
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

sns.heatmap(
    covmat_new / np.outer(np.sqrt(np.diag(covmat_new)), np.sqrt(np.diag(covmat_new))),
    annot=False,
    cmap="YlGnBu",
    ax=axs[0],
    cbar_ax=cbar_ax,
)

sns.heatmap(
    covmat_old / np.outer(np.sqrt(np.diag(covmat_old)), np.sqrt(np.diag(covmat_old))),
    annot=False,
    cmap="YlGnBu",
    ax=axs[1],
    cbar_ax=cbar_ax,
)

axs[0].set_title("Legacy Correlation Matrix")
axs[1].set_title("New Correlation Matrix")
plt.show()