"""
Test to ensure the validphys.pseudodata.get_pseudodata action
correctly obtains the appropriate pseudodata for an n3fit fit.

A fit has been generated called pseudodata_test_fit_n3fit
which has the pseudodata saved as training and validation splits.
This is used to benchmark the correctness of the pseudodata
recreation.
"""
import pandas as pd
import numpy as np
import pytest

from validphys.api import API
from validphys.tests.conftest import FIT, PSEUDODATA_FIT


def test_read_fit_pseudodata():
    fit_pseudodata = API.read_fit_pseudodata(fit=PSEUDODATA_FIT)

    nrep = API.num_fitted_replicas(fit=PSEUDODATA_FIT)
    assert nrep == len(fit_pseudodata)

    for data, tr_idx, val_idx in fit_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_read_pdf_pseudodata():
    pdf_pseudodata = API.read_pdf_pseudodata(fit=PSEUDODATA_FIT)

    pdf = API.pdf(pdf=PSEUDODATA_FIT)
    # -1 because we ignore replica 0
    assert len(pdf) - 1 == len(pdf_pseudodata)

    for data, tr_idx, val_idx in pdf_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_recreate_fit_pseudodata():
    fit_pseudodata = API.recreate_fit_pseudodata(fit=PSEUDODATA_FIT)

    nrep = API.num_fitted_replicas(fit=PSEUDODATA_FIT)
    assert nrep == len(fit_pseudodata)

    for data, tr_idx, val_idx in fit_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_recreate_pdf_pseudodata():
    pdf_pseudodata = API.recreate_pdf_pseudodata(fit=PSEUDODATA_FIT)

    pdf = API.pdf(pdf=PSEUDODATA_FIT)
    # -1 because we ignore replica 0
    assert len(pdf) - 1 == len(pdf_pseudodata)

    for data, tr_idx, val_idx in pdf_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_no_savepseudodata():
    for func in (API.read_fit_pseudodata, API.read_pdf_pseudodata):
        with pytest.raises(FileNotFoundError):
            # Check a FileNotFoundError is raised
            # if the input fit wasn't generated
            # with the savepseudodata flag set to true
            func(fit=FIT)


def test_read_matches_recreate():
    reads = API.read_fit_pseudodata(fit=PSEUDODATA_FIT)
    recreates = API.recreate_fit_pseudodata(fit=PSEUDODATA_FIT)
    for read, recreate in zip(reads, recreates):
        # We ignore the absolute ordering of the dataframes and just check
        # that they contain identical elements.
        pd.testing.assert_frame_equal(
            read.pseudodata, recreate.pseudodata, check_like=True
        )
        pd.testing.assert_index_equal(read.tr_idx, recreate.tr_idx, check_order=False)
        pd.testing.assert_index_equal(read.val_idx, recreate.val_idx, check_order=False)


def test_make_level0_data():
    
    dataset='NMC'
    pdfname='NNPDF40_nnlo_as_01180'
    theoryid=200

    l0_cd = API.make_level0_data(dataset_inputs = [{"dataset":dataset}], 
                                use_cuts="internal", theoryid=theoryid, fakepdf = pdfname)
                                
    l0_vals = l0_cd[0].central_values
    t0_pred = np.array([0.3545763 , 0.36826616, 0.36873821, 0.37127424, 0.37290897,
       0.37316775, 0.37055304, 0.3722504 , 0.37302464, 0.36924048,
       0.37034822, 0.37095343, 0.36640187, 0.36637822, 0.36574894,
       0.35825978, 0.35653348, 0.35473318, 0.35365563, 0.34134523,
       0.33823255, 0.33386657, 0.33338396, 0.31109264, 0.30324802,
       0.30229506, 0.25416182, 0.2523507 , 0.37196013, 0.37330645,
       0.37654302, 0.37939536, 0.38006708, 0.36966435, 0.37476681,
       0.37801795, 0.38023252, 0.38028246, 0.36852294, 0.37221565,
       0.37471784, 0.37672207, 0.37771963, 0.37009838, 0.37165855,
       0.37282795, 0.37345963, 0.37297594, 0.36651327, 0.36682826,
       0.36658552, 0.36602139, 0.36449717, 0.35922251, 0.35652948,
       0.35558624, 0.3538604 , 0.35223064, 0.34976281, 0.33977084,
       0.33597368, 0.33219177, 0.33006355, 0.32609279, 0.30810523,
       0.29984762, 0.2973673 , 0.29337914, 0.28835707, 0.26067163,
       0.2463048 , 0.23991642, 0.2346202 , 0.15752154, 0.14484311,
       0.37357085, 0.37989003, 0.38811387, 0.39269338, 0.39623239,
       0.38779016, 0.39516254, 0.39979762, 0.40095537, 0.37504614,
       0.38304588, 0.38946276, 0.3957882 , 0.39952722, 0.39975423,
       0.37745691, 0.38209647, 0.3865788 , 0.39086245, 0.39332198,
       0.39323533, 0.37680676, 0.37996398, 0.38281809, 0.38475522,
       0.38561046, 0.38429377, 0.37467448, 0.37630849, 0.37726671,
       0.37769016, 0.37684543, 0.36781676, 0.36746375, 0.36715733,
       0.36632717, 0.36472775, 0.36155728, 0.35540615, 0.35320503,
       0.35029904, 0.34840341, 0.34560717, 0.34203748, 0.33253254,
       0.32686904, 0.32281255, 0.31914235, 0.31415925, 0.29386229,
       0.29031662, 0.28439739, 0.27832586, 0.27371331, 0.24110857,
       0.23133557, 0.22451033, 0.21709159, 0.14911905, 0.13010754,
       0.12310468, 0.39012269, 0.39390918, 0.39924528, 0.40293211,
       0.4056547 , 0.39786495, 0.40730062, 0.41398714, 0.41582426,
       0.41570175, 0.39194794, 0.40105368, 0.40975066, 0.41707675,
       0.41820279, 0.41567735, 0.39337477, 0.39929797, 0.40673509,
       0.4117446 , 0.41284718, 0.40930253, 0.39081395, 0.3929952 ,
       0.39747039, 0.40079244, 0.40128669, 0.39749262, 0.38580429,
       0.38676072, 0.38910638, 0.39044752, 0.38910755, 0.37872247,
       0.37961559, 0.38016427, 0.3793288 , 0.37651091, 0.36785649,
       0.36804973, 0.36655191, 0.36481613, 0.3620474 , 0.34931946,
       0.3463118 , 0.34349012, 0.34000501, 0.33611505, 0.3238292 ,
       0.31964367, 0.31530522, 0.31121644, 0.30637345, 0.29586067,
       0.28526531, 0.27994941, 0.27426498, 0.26944325, 0.26316807,
       0.24346269, 0.22612249, 0.21976363, 0.21375908, 0.20546207,
       0.15822622, 0.13056665, 0.12244491, 0.11192735])
    

    assert(np.abs(np.sum(t0_pred / l0_vals) / len(l0_vals) - 1) <= 1e-9)