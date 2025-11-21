"""
test_loader.py

Test overfit metric implementation.
"""

from validphys.api import API
from validphys.tests.conftest import FIT_3REPLICAS, FIT_3REPLICAS_DIAG

config = {
    "use_t0": True,
    "use_cuts": "fromfit",
    "theory": {"from_": "fit"},
    "theoryid": {"from_": "theory"},
    "datacuts": {"from_": "fit"},
    "t0pdfset": {"from_": "datacuts"},
    "pdf": {"from_": "fit"},
    "dataset_inputs": {"from_": "fit"},
}


def test_overfit_chi2():
    for fit in [FIT_3REPLICAS]:
        # TODO: implement diagonal basis test here. For this we should save the rotation matrix
        #  as part of the pseudodata to reconstruct the vl chi2
        diagonal_basis = True if fit == FIT_3REPLICAS_DIAG else False
        replica_info = API.replica_data(**config, fit=fit, diagonal_basis=diagonal_basis)
        val_erf = [info.validation for info in replica_info]
        chi2s_per_replica = API.calculate_chi2s_per_replica(
            **config, fit=fit, diagonal_basis=diagonal_basis
        )
        assert (abs(chi2s_per_replica.diagonal() - val_erf) < 1e-4).all()
