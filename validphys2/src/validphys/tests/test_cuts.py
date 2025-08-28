from validphys.api import API
from validphys.core import InternalCutsWrapper, SimilarCuts
from validphys.tests.conftest import DATA, PDF, THEORYID


def test_similarity_cuts():
    plain = [{"dataset": dt["dataset"]} for dt in DATA]
    inp = {
        "theoryid": THEORYID,
        "pdf": PDF,
        "cut_similarity_threshold": 1.5,
        "use_cuts": "fromsimilarpredictions",
        "cuts_intersection_spec": [{"dataset_inputs": DATA}, {"dataset_inputs": plain}],
        "dataset_input": DATA[1],
    }
    ds = API.dataset(**inp)
    assert isinstance(ds.cuts, SimilarCuts)

    inp["do_not_require_similarity_for"] = [DATA[1]["dataset"]]
    ds = API.dataset(**inp)
    assert isinstance(ds.cuts, InternalCutsWrapper)
