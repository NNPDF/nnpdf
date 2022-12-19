from validphys.loader import Loader
from validphys.tests.conftest import SINGLE_DATASET
from validphys.inconsys import process_commondata

def test_process_commondata():
    l = Loader()
    obs = SINGLE_DATASET["dataset"]
    commondata = l.check_commondata(obs).load_commondata_instance()
    # test that if dataset is not within inconsistent_datasets, nothing is changed
    commondata_1 = process_commondata(commondata=commondata,ADD=False,MULT=False,
                        CORR=False,UNCORR=False,inconsistent_datasets=[],sys_rescaling_factor=2)

    assert( commondata == commondata_1)

    # check that we create a new commondata instance when dataset is within inconsistent_datasets
    commondata_2 = process_commondata(commondata=commondata,ADD=True,MULT=True,
                        CORR=False,UNCORR=False,inconsistent_datasets=[obs],sys_rescaling_factor=1)

    assert( commondata != commondata_2)