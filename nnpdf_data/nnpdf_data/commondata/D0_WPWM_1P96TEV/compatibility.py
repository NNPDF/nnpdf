from validphys.api import API
import numpy as np

new_implementation = "D0_WPWM_1P96TEV_ASY"
old_implementation = "D0_WPWM_1P96TEV_ASY"

inp1 = {"dataset_input": {"dataset": f"{new_implementation}"}, "theoryid": 40_000_000, "use_cuts": "internal", "t0pdfset": "NNPDF40_nnlo_as_01180", "use_t0": True}
inp2 = {"dataset_input": {"dataset": f"{old_implementation}", "variant": "legacy"}, "theoryid": 40_000_000, "use_cuts": "internal", "t0pdfset": "NNPDF40_nnlo_as_01180", "use_t0": True}

covmat1 = API.covmat_from_systematics(**inp1)
covmat2 = API.covmat_from_systematics(**inp2)

t0_covmat1 = API.t0_covmat_from_systematics(**inp1)
t0_covmat2 = API.t0_covmat_from_systematics(**inp2)

print(np.argwhere(~np.isclose(covmat1, covmat2)))
print(np.argwhere(~np.isclose(t0_covmat1, t0_covmat2)))