"""
    Library of function for reading  NNPDF objects
"""
import hashlib
import copy
import numpy as np

from NNPDF import RandomGenerator
from validphys.core import ExperimentSpec as vp_Exp
from validphys.core import DataSetSpec as vp_Dataset


def make_tr_val_mask(datasets, exp_name, seed):
    """
    Generates the training and validation masks for a given experiment

    # Arguments:
        - `datasets`: list of datasets corresponding to one experiment
        - `exp_name`: name of the experiment
        - `seed` : seed for the random tr/vl split

    # Return
        - `trmask`: numpy boolean array
        - `vlmask`: not trmask
    """
    # Set the seed for the experiment
    nameseed = int(hashlib.sha256(exp_name.encode()).hexdigest(), 16) % 10 ** 8
    nameseed += seed
    np.random.seed(nameseed)

    trmask_partial = []
    vlmask_partial = []
    for dataset_dict in datasets:
        ndata = dataset_dict["ndata"]
        frac = dataset_dict["frac"]
        trmax = int(frac * ndata)
        mask = np.concatenate(
            [np.ones(trmax, dtype=np.bool), np.zeros(ndata - trmax, dtype=np.bool)]
        )
        np.random.shuffle(mask)
        trmask_partial.append(mask)
        vlmask_partial.append(mask == False)
    trmask = np.concatenate(trmask_partial)
    vlmask = np.concatenate(vlmask_partial)

    return trmask, vlmask


def fk_parser(fk, is_hadronic=False):
    """
    # Arguments:
        - `fk`: fktable object

    # Return:
        - `dict_out`: dictionary with all information about the fktable
            - 'xgrid'
            - 'nx'
            - 'ndata'
            - 'basis'
            - 'fktable'
    """

    # Dimensional information
    nonzero = fk.GetNonZero()
    ndata = fk.GetNData()
    nx = fk.GetTx()

    # Basis of active flavours
    nbasis = nonzero
    basis = fk.get_flmap()

    # Read the xgrid and the fktable to numpy arrays
    xgrid_flat = fk.get_xgrid()
    fktable_flat = fk.get_sigma()

    # target shape
    if is_hadronic:
        nxsqrt = int(np.sqrt(nx))
        shape_out = (ndata, nbasis, nxsqrt, nxsqrt)
        xgrid = xgrid_flat.reshape(1, nxsqrt)
    else:
        shape_out = (ndata, nbasis, nx)
        xgrid = xgrid_flat.reshape(1, nx)

    # remove padding from the fktable (if any)
    l = len(fktable_flat)
    # get the padding
    pad_position = fk.GetDSz()
    pad_size = fk.GetDSz() - nx * nonzero
    # remove the padding
    if pad_size > 0:
        mask = np.ones(l, dtype=bool)
        for i in range(1, int(l / pad_position) + 1):
            marker = pad_position * i
            mask[marker - pad_size : marker] = False
        fktable_array = fktable_flat[mask]
    else:
        fktable_array = fktable_flat
    # reshape
    fktable = fktable_array.reshape(shape_out)

    dict_out = {
        "ndata": ndata,
        "nbasis": nbasis,
        "nonzero": nbasis,
        "basis": basis,
        "nx": nx,
        "xgrid": xgrid,
        "fktable": fktable,
    }
    return dict_out


def common_data_reader_dataset(dataset_c, dataset_spec):
    """
    Import fktable, common data and experimental data for the given data_name

    # Arguments:
        - `dataset_c`: c representation of the dataset object
        - `dataset_spec`: python representation of the dataset object

    #Returns:
        - `[dataset_dict]`: a (len 1 list of) dictionary with:
            - dataset: full dataset object from which all data has been extracted
            - hadronic: is hadronic?
            - operation: operation to perform on the fktables below
            - fktables : a list of dictionaries with the following items
                - ndata: number of data points
                - nonzero/nbasis: number of PDF entries which are non zero
                - basis: array containing the information of the non-zero PDF entries
                - nx: number of points in the xgrid (1-D), i.e., for hadronic the total number is nx * nx
                - xgrid: grid of x points
                - fktable: 3/4-D array of shape (ndata, nonzero, nx, (nx))

    instead of the dictionary object that model_gen needs
    """
    how_many = dataset_c.GetNSigma()
    dict_fktables = []
    for i in range(how_many):
        fktable = dataset_c.GetFK(i)
        dict_fktables.append(fk_parser(fktable, dataset_c.IsHadronic()))

    dataset_dict = {
        "fktables": dict_fktables,
        "hadronic": dataset_c.IsHadronic(),
        "operation": dataset_spec.op,
        "name": dataset_c.GetSetName(),
        "frac": dataset_spec.frac,
        "ndata": dataset_c.GetNData(),
    }

    return [dataset_dict]


def common_data_reader_experiment(experiment_c, experiment_spec):
    """
        Wrapper around the experiments. Loop over all datasets in an experiment,
        calls common_data_reader on them and return a list with the content.

        # Arguments:
            - `experiment_c`: c representation of the experiment object
            - `experiment_spec`: python representation of the experiment object

        # Returns:
            - `[parsed_datasets]`: a list of dictionaries output from `common_data_reader_dataset`
    """
    parsed_datasets = []
    for dataset_c, dataset_spec in zip(
        experiment_c.DataSets(), experiment_spec.datasets
    ):
        parsed_datasets += common_data_reader_dataset(dataset_c, dataset_spec)
    return parsed_datasets


def common_data_reader(
    spec,
    t0pdfset,
    replica_seeds=None,
    trval_seeds=None,
    kpartitions=None,
    rotate_diagonal=False
):
    """
    Wrapper to read the information from validphys object
    This function receives either a validphyis experiment or dataset objects

    # Returns:
        - `all_dict_out`: a dictionary containing all the information of the experiment/dataset
        for training, validation and experimental
                'datasets' : list of the datasets contained in the experiment
                'name' : name of the experiment
                'expdata_true' : non-replica data
                'invcovmat_true' : inverse of the covmat (non-replica)

                'trmask' : mask for the training data
                'invcovmat' : inverse of the covmat for the training data
                'ndata' : number of datapoints for the training data
                'expdata' : experimental data (replica'd) for training

                'vlmask' : (same as above for validation)
                'invcovmat_vl' : (same as above for validation)
                'ndata_vl' : (same as above for validation)
                'expdata_vl' :  (same as above for validation)

                'positivity' : bool - is this a positivity set?
                'count_chi2' : should this be counted towards the chi2
    """
    if replica_seeds is None:
        replica_seeds = []
    if trval_seeds is None:
        trval_seeds = [0]
    if kpartitions is None:
        kpartitions = []
    # TODO
    # This whole thing would be much more clear / streamlined if
    #   - The c experiment/dataset object had all the required information for the fit
    #       (i.e., there is a swig conversion for everything, right now missing the operator)
    #   - The python object stored the load within the spec when doing spec.load()
    #                                   this way it would not be necessary to load twice
    #   - The python object had all necessary information (same as point 1 but inverted)

    spec_c = spec.load()
    ndata = spec_c.GetNData()
    expdata_true = spec_c.get_cv().reshape(1, ndata)
    if t0pdfset is not None:
        spec_c.SetT0(t0pdfset)
    base_mcseed = int(hashlib.sha256(str(spec).encode()).hexdigest(), 16) % 10 ** 8

    if replica_seeds:
        all_expdatas = []
    else:
        all_expdatas = [expdata_true.reshape(ndata)]

    for replica_seed in replica_seeds:
        spec_replica_c = type(spec_c)(spec_c) # I might need the t0 set here as well

        # Replica generation
        mcseed = base_mcseed + replica_seed
        RandomGenerator.InitRNG(0, mcseed)
        spec_replica_c.MakeReplica()
        all_expdatas.append(spec_replica_c.get_cv())

    # spec_c = spec_replica_c

    if isinstance(spec, vp_Exp):
        datasets = common_data_reader_experiment(spec_c, spec)
    elif isinstance(spec, vp_Dataset):
        datasets = common_data_reader_dataset(spec_c, spec)
    else:
        raise ValueError(
            "reader.py: common_data_reader, didn't understood spec type: {0}".format(
                type(spec)
            )
        )

    folds = {"training": [], "validation": [], "experimental" : []}
    for partition in kpartitions:
        data_fold = partition.get("datasets", [])
        mask = []
        for dataset in datasets:
            if dataset['name'] in data_fold:
                mask.append(np.zeros(dataset['ndata'], dtype=np.bool))
            else:
                mask.append(np.ones(dataset['ndata'], dtype=np.bool))
        folds["experimental"].append(np.concatenate(mask))

    exp_name = spec.name
    covmat = spec_c.get_covmat()
    inv_true = np.linalg.inv(covmat)

    if rotate_diagonal:
        eig, v = np.linalg.eigh(covmat)
        dt_trans = v.T
    else:
        dt_trans = None

    # Now it is time to build the masks for the training validation split
    all_dict_out = []
    for expdata, trval_seed in zip(all_expdatas, trval_seeds):
        tr_mask, vl_mask = make_tr_val_mask(datasets, exp_name, seed=trval_seed)

        if rotate_diagonal:
            expdata = np.matmul(dt_trans, expdata)
            # make a 1d array of the diagonal
            covmat_tr = eig[tr_mask]
            invcovmat_tr = 1./covmat_tr

            covmat_vl = eig[vl_mask]
            invcovmat_vl = 1./covmat_vl
        else:
            covmat_tr = covmat[tr_mask].T[tr_mask]
            invcovmat_tr = np.linalg.inv(covmat_tr)

            covmat_vl = covmat[vl_mask].T[vl_mask]
            invcovmat_vl = np.linalg.inv(covmat_vl)

        ndata_tr = np.count_nonzero(tr_mask)
        expdata_tr = expdata[tr_mask].reshape(1, ndata_tr)

        ndata_vl = np.count_nonzero(vl_mask)
        expdata_vl = expdata[vl_mask].reshape(1, ndata_vl)


        for ex_fold in folds["experimental"]:
            folds["training"].append(ex_fold[tr_mask])
            folds["validation"].append(ex_fold[vl_mask])

        dict_out = {
            "datasets": datasets,
            "name": exp_name,
            "expdata_true": expdata_true,
            "invcovmat_true": inv_true,
            "trmask": tr_mask,
            "invcovmat": invcovmat_tr,
            "ndata": ndata_tr,
            "expdata": expdata_tr,
            "vlmask": vl_mask,
            "invcovmat_vl": invcovmat_vl,
            "ndata_vl": ndata_vl,
            "expdata_vl": expdata_vl,
            "positivity": False,
            "count_chi2": True,
            "folds" : folds,
            "data_transformation": dt_trans,
        }
        all_dict_out.append(dict_out)

    return all_dict_out


def positivity_reader(pos_spec):
    """
    Specific reader for positivity sets
    """
    pos_c = pos_spec.load()
    ndata = pos_c.GetNData()

    parsed_set = [fk_parser(pos_c, pos_c.IsHadronic())]

    pos_sets = [
        {
            "fktables": parsed_set,
            "hadronic": pos_c.IsHadronic(),
            "operation": "NULL",
            "name": pos_spec.name,
            "frac": 1.0,
            "ndata": ndata,
        }
    ]

    positivity_factor = pos_spec.poslambda

    dict_out = {
        "datasets": pos_sets,
        "trmask": np.ones(ndata, dtype=np.bool),
        "name": pos_spec.name,
        "expdata": np.zeros((1, ndata)),
        "ndata": ndata,
        "positivity": True,
        "lambda": positivity_factor,
        "count_chi2": False,
    }

    return dict_out
