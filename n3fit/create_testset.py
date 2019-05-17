import numpy as np
import os

from reportengine.compat import yaml

def test_runcard(runcard_file, datasets_out):
    """
        Given a runcard and a set of datasets, generate a TEST-runcard where the datasets in datasets_out are removed
        from the training and put together as a new experiment called TEST
    """
    # Load the input runcard as a dictionary
    f = open(runcard_file, 'r')
    runcard_dict = yaml.safe_load(f)
    runcard_exp = runcard_dict['experiments']
    f.close()

    # Generate the list of datasets which are to be left out
    data_tests = []
    for dataset in datasets_out:
        data_tests.append(
                {'dataset' : dataset, 'frac' : 1.0} # TODO copy the dictionary from the runcard instead
                )

    # Create a test experiment
    from ModelTrainer import TESTNAME
    test_experiment = {
            'experiment' : TESTNAME,
            'datasets' : data_tests
            }

    # Now drop experiments and datasets that are being used for training
    for experiment in runcard_exp:
        datasets = list(filter(lambda x: x['dataset'] not in datasets_out, experiment['datasets']))
        if datasets:
            experiment['datasets'] = datasets
        else:
            runcard_exp.remove(experiment)

    # Generate the new runcard 
    # TODO: maybe there is a proper report-enginy way of doing this
    runcard_exp.append(test_experiment)
    new_runcard = "TEST-{0}".format(os.path.basename(runcard_file))
    no = open(new_runcard, 'w')
    yaml.round_trip_dump(runcard_dict, no, explicit_start=True, explicit_end=True,
                  default_flow_style=True)
    no.close()

    print("\n > You can find your new runcard at {0}/{1}\n".format(os.getcwd(), new_runcard))

def create_testset(experiments, runcard_file = "runcards/NNPDF31_nlo_as_0118.yml"):
    """ 
        This function gets the list of experiment and goes through the datasets
        reading their experiment type and range in x
        Depending on these, chooses the ones to leave out for testing
    """

    # 1) Read all datasets in the runcard and classify them by type
    dataset_by_proc = {}
    all_experiments = {}
    for exp in experiments:
        all_experiments[str(exp)] = len(exp.datasets)
        for dataset in exp.datasets:
            data_c = dataset.load()
            proc_type = data_c.GetProc(0)
            fktable = data_c.GetFK(0)
            xgrid = fktable.get_xgrid()
            xmin = np.min(xgrid)
            xmax = np.max(xgrid)
            data_dict = {
                    'name' : str(dataset),
                    'xmin' : xmin,
                    'xmax' : xmax,
                    'exp' : str(exp),
                    }
            if proc_type in dataset_by_proc.keys():
                dataset_by_proc[proc_type].append(data_dict)
            else:
                dataset_by_proc[proc_type] = [data_dict]

    select_min = False # If true selects the datasets with the smallest min(x), recommended: False

    # 2) Now for every process type with more than one dataset, leave out the smallest
    datasets_out = []
    for proc_type, datasets in dataset_by_proc.items():
        l = len(datasets)
        if l != 1:
            xmin = 1.0
            xmax = 0.0
            for dataset in datasets:
                xm = dataset['xmin']
                if select_min:
                    if xm < xmin:
                        xmin = xmin
                        data_out = dataset
                else:
                    if xm > xmax:
                        xmax = xm
                        data_out = dataset
            print("For {0} we find {1} datasets, leaving one out: {2}".format(proc_type, l, data_out))
            datasets_out.append(data_out['name'])
            all_experiments[dataset['exp']] -= 1

    if len(datasets_out) == 0:
        print("No redundant datasets were found")
    else:
        test_runcard(runcard_file, datasets_out)
