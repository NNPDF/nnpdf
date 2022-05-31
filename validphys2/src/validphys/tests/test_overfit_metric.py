"""
test_loader.py

Test overfit metric implementation.
"""

from validphys.api import API
from validphys.tests.conftest import FIT_3REPLICAS

config = {
    "fit": FIT_3REPLICAS,
    "use_t0": True,
    "use_cuts": "fromfit",
    "theory": {"from_": "fit"},
    "theoryid": {"from_": "theory"},
    "datacuts": {"from_": "fit"},
    "t0pdfset": {"from_": "datacuts"},
    "pdfs": [{"from_": "fit"}],
    "dataset_inputs": {"from_": "fit"},
}


def test_overfit_chi2():
    replica_info = API.replica_data(**config)
    val_erf = [info.validation for info in replica_info]
    chi2s_per_replica = API.calculate_chi2s_per_replica(**config)
    assert (abs(chi2s_per_replica.diagonal() - val_erf) < 1e-4).all()
