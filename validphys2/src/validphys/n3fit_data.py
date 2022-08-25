"""
n3fit_data.py

Providers which prepare the data ready for
:py:func:`n3fit.performfit.performfit`. Returns python objects but the underlying
functions make calls to libnnpdf C++ library.

"""
import functools
from collections import defaultdict
import hashlib
import logging

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table

from validphys.n3fit_data_utils import (
    validphys_group_extractor,
)
from validphys.core import IntegrabilitySetSpec, TupleComp

log = logging.getLogger(__name__)


def replica_trvlseed(replica, trvlseed, same_trvl_per_replica=False):
    """Generates the ``trvlseed`` for a ``replica``."""
    # TODO: move to the new infrastructure
    # https://numpy.org/doc/stable/reference/random/index.html#introduction
    np.random.seed(seed=trvlseed)
    if same_trvl_per_replica:
        return np.random.randint(0, pow(2, 31))
    for _ in range(replica):
        res = np.random.randint(0, pow(2, 31))
    return res


def replica_nnseed(replica, nnseed):
    """Generates the ``nnseed`` for a ``replica``."""
    np.random.seed(seed=nnseed)
    for _ in range(replica):
        res = np.random.randint(0, pow(2, 31))
    return res


def replica_mcseed(replica, mcseed, genrep):
    """Generates the ``mcseed`` for a ``replica``."""
    if not genrep:
        return None
    np.random.seed(seed=mcseed)
    for _ in range(replica):
        res = np.random.randint(0, pow(2, 31))
    return res


class _TrMasks(TupleComp):
    """Class holding the training validation mask for a group of datasets
    If the same group of dataset receives the same trvlseed then the mask
    will be the same.
    This class holds said information so it can be reused easily, i.e.,
    ``group_name`` and ``seed`` define the ``masks``.
    """

    def __init__(self, group_name, seed, masks=None):
        self.masks = masks
        super().__init__(group_name, seed)

    def __iter__(self):
        for m in self.masks:
            yield m


def tr_masks(data, replica_trvlseed):
    """Generate the boolean masks used to split data into training and
    validation points. Returns a list of 1-D boolean arrays, one for each
    dataset. Each array has length equal to N_data, the datapoints which
    will be included in the training are ``True`` such that

        tr_data = data[tr_mask]

    """
    nameseed = int(hashlib.sha256(str(data).encode()).hexdigest(), 16) % 10**8
    nameseed += replica_trvlseed
    # TODO: update this to new random infrastructure.
    np.random.seed(nameseed)
    trmask_partial = []
    for dataset in data.datasets:
        # TODO: python commondata will not require this rubbish.
        # all data if cuts are None
        cuts = dataset.cuts
        ndata = len(cuts.load()) if cuts else dataset.commondata.ndata
        frac = dataset.frac
        trmax = int(frac * ndata)
        mask = np.concatenate(
            [np.ones(trmax, dtype=np.bool), np.zeros(ndata - trmax, dtype=np.bool)]
        )
        np.random.shuffle(mask)
        trmask_partial.append(mask)
    return _TrMasks(str(data), replica_trvlseed, trmask_partial)


def kfold_masks(kpartitions, data):
    """Collect the masks (if any) due to kfolding for this data.
    These will be applied to the experimental data before starting
    the training of each fold.

    Parameters
    ----------
    kpartitions: list[dict]
        list of partitions, each partition dictionary with key-value pair
        `datasets` and a list containing the names of all datasets in that
        partition. See n3fit/runcards/Basic_hyperopt.yml for an example
        runcard or the hyperopt documentation for an expanded discussion on
        k-fold partitions.
    data: validphys.core.DataGroupSpec
        full list of data which is to be partitioned.

    Returns
    -------
    kfold_masks: list[np.array]
        A list containing a boolean array for each partition. Each array is
        a 1-D boolean array with length equal to the number of cut datapoints
        in ``data``. If a dataset is included in a particular fold then the
        mask will be True for the elements corresponding to those datasets
        such that data.load().get_cv()[kfold_masks[i]] will return the
        datapoints in the ith partition. See example below.

    Examples
    --------
    >>> from validphys.api import API
    >>> partitions=[
    ...     {"datasets": ["HERACOMBCCEM", "HERACOMBNCEP460", "NMC", "NTVNBDMNFe"]},
    ...     {"datasets": ["HERACOMBCCEP", "HERACOMBNCEP575", "NMCPD", "NTVNUDMNFe"]}
    ... ]
    >>> ds_inputs = [{"dataset": ds} for part in partitions for ds in part["datasets"]]
    >>> kfold_masks = API.kfold_masks(dataset_inputs=ds_inputs, kpartitions=partitions, theoryid=53, use_cuts="nocuts")
    >>> len(kfold_masks) # one element for each partition
    2
    >>> kfold_masks[0] # mask which splits data into first partition
    array([False, False, False, ...,  True,  True,  True])
    >>> data = API.data(dataset_inputs=ds_inputs, theoryid=53, use_cuts="nocuts")
    >>> fold_data = data.load().get_cv()[kfold_masks[0]]
    >>> len(fold_data)
    604
    >>> kfold_masks[0].sum()
    604

    """
    list_folds = []
    if kpartitions is not None:
        for partition in kpartitions:
            data_fold = partition.get("datasets", [])
            mask = []
            for dataset in data.datasets:
                # TODO: python commondata will not require this rubbish.
                # all data if cuts are None
                cuts = dataset.cuts
                ndata = len(cuts.load()) if cuts else dataset.commondata.ndata
                # If the dataset is in the fold, its mask is full of 0s
                if str(dataset) in data_fold:
                    mask.append(np.zeros(ndata, dtype=np.bool))
                # otherwise of ones
                else:
                    mask.append(np.ones(ndata, dtype=np.bool))
            list_folds.append(np.concatenate(mask))
    return list_folds


@functools.lru_cache
def fittable_datasets_masked(data, tr_masks):
    """Generate a list of :py:class:`validphys.n3fit_data_utils.FittableDataSet`
    from a group of dataset and the corresponding training/validation masks
    """
    # This is separated from fitting_data_dict so that we can cache the result
    # when the trvlseed is the same for all replicas (great for parallel replicas)
    return validphys_group_extractor(data.datasets, tr_masks.masks)


def fitting_data_dict(
    data,
    make_replica,
    dataset_inputs_loaded_cd_with_cuts,
    dataset_inputs_fitting_covmat,
    tr_masks,
    kfold_masks,
    fittable_datasets_masked,
    diagonal_basis=None,
):
    """
    Provider which takes  the information from validphys ``data``.

    Returns
    -------
    all_dict_out: dict
        Containing all the information of the experiment/dataset
        for training, validation and experimental With the following keys:

        'datasets'
            list of dictionaries for each of the datasets contained in ``data``
        'name'
            name of the ``data`` - typically experiment/group name
        'expdata_true'
            non-replica data
        'covmat'
            full covmat
        'invcovmat_true'
            inverse of the covmat (non-replica)
        'trmask'
            mask for the training data
        'invcovmat'
            inverse of the covmat for the training data
        'ndata'
            number of datapoints for the training data
        'expdata'
            experimental data (replica'd) for training
        'vlmask'
            (same as above for validation)
        'invcovmat_vl'
            (same as above for validation)
        'ndata_vl'
            (same as above for validation)
        'expdata_vl'
            (same as above for validation)
        'positivity'
            bool - is this a positivity set?
        'count_chi2'
            should this be counted towards the chi2
    """
    # TODO: Plug in the python data loading when available. Including but not
    # limited to: central values, ndata, replica generation, covmat construction
    expdata_true = np.concatenate([d.central_values for d in dataset_inputs_loaded_cd_with_cuts])

    expdata = make_replica
    tr_masks = tr_masks.masks
    covmat = dataset_inputs_fitting_covmat # t0 covmat, or theory covmat or whatever was decided by the runcard
    inv_true = np.linalg.inv(covmat)
    fittable_datasets = fittable_datasets_masked

    if diagonal_basis:
        log.info("working in diagonal basis.")
        eig, v = np.linalg.eigh(covmat)
        dt_trans = v.T
    else:
        dt_trans = None
        dt_trans_tr = None
        dt_trans_vl = None

    tr_mask = np.concatenate(tr_masks)
    vl_mask = ~tr_mask

    if diagonal_basis:
        expdata = np.matmul(dt_trans, expdata)
        # make a 1d array of the diagonal
        covmat_tr = eig[tr_mask]
        invcovmat_tr = 1.0 / covmat_tr

        covmat_vl = eig[vl_mask]
        invcovmat_vl = 1.0 / covmat_vl

        # prepare a masking rotation
        dt_trans_tr = dt_trans[tr_mask]
        dt_trans_vl = dt_trans[vl_mask]
    else:
        covmat_tr = covmat[tr_mask].T[tr_mask]
        invcovmat_tr = np.linalg.inv(covmat_tr)

        covmat_vl = covmat[vl_mask].T[vl_mask]
        invcovmat_vl = np.linalg.inv(covmat_vl)

    ndata_tr = np.count_nonzero(tr_mask)
    expdata_tr = expdata[tr_mask].reshape(1, ndata_tr)

    ndata_vl = np.count_nonzero(vl_mask)
    expdata_vl = expdata[vl_mask].reshape(1, ndata_vl)

    # Now save a dictionary of training/validation/experimental folds
    # for training and validation we need to apply the tr/vl masks
    # for experimental we need to negate the mask
    folds = defaultdict(list)
    for fold in kfold_masks:
        folds["training"].append(fold[tr_mask])
        folds["validation"].append(fold[vl_mask])
        folds["experimental"].append(~fold)

    # This dictionary contains a list of fittable datasets
    # which contains the instructions on how to generate each observable for the fit
    # plus the information that glue all of them together (covmat, ndata, etc)
    # TODO: for consistency with the rest of validphys a FittableGroup should be created
    dict_out = {
        "datasets": fittable_datasets,
        "name": str(data),
        "expdata_true": expdata_true.reshape(1, -1),
        "invcovmat_true": inv_true,
        "covmat": covmat,
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
        "folds": folds,
        "data_transformation_tr": dt_trans_tr,
        "data_transformation_vl": dt_trans_vl,
    }
    return dict_out


exps_fitting_data_dict = collect("fitting_data_dict", ("group_dataset_inputs_by_metadata",))


def replica_nnseed_fitting_data_dict(replica, exps_fitting_data_dict, replica_nnseed):
    """For a single replica return a tuple of the inputs to this function.
    Used with `collect` over replicas to avoid having to perform multiple
    collects.

    See Also
    --------
    replicas_nnseed_fitting_data_dict - the result of collecting this function
    over replicas.

    """
    return (replica, exps_fitting_data_dict, replica_nnseed)


replicas_nnseed_fitting_data_dict = collect("replica_nnseed_fitting_data_dict", ("replicas",))
groups_replicas_indexed_make_replica = collect(
    "indexed_make_replica", ("group_dataset_inputs_by_experiment", "replicas")
)


@table
def pseudodata_table(groups_replicas_indexed_make_replica, replicas):
    """Creates a pandas DataFrame containing the generated pseudodata. The
    index is :py:func:`validphys.results.experiments_index` and the columns
    are the replica numbers.

    Notes
    -----
    Whilst running ``n3fit``, this action will only be called if
    `fitting::savepseudodata` is `true` and replicas are fitted one at a time.
    The table can be found in the replica folder i.e. <fit dir>/nnfit/replica_*/

    """
    # Concatenate over replicas
    df = pd.concat(groups_replicas_indexed_make_replica)
    df.columns = [f"replica {rep}" for rep in replicas]
    return df


@table
def training_pseudodata(pseudodata_table, training_mask):
    """Save the training data for the given replica.
    Activate by setting ``fitting::savepseudodata: True``
    from within the fit runcard.

    See Also
    --------
    :py:func:`validphys.n3fit_data.validation_pseudodata`
    """
    return pseudodata_table.loc[training_mask.values]


@table
def validation_pseudodata(pseudodata_table, training_mask):
    """Save the training data for the given replica.
    Activate by setting ``fitting::savepseudodata: True``
    from within the fit runcard.

    See Also
    --------
    :py:func:`validphys.n3fit_data.training_pseudodata`
    """
    return pseudodata_table.loc[~training_mask.values]


exps_tr_masks = collect("tr_masks", ("group_dataset_inputs_by_experiment",))
replicas_exps_tr_masks = collect("exps_tr_masks", ("replicas",))


@table
def replica_training_mask_table(replica_training_mask):
    """Same as ``replica_training_mask`` but with a table decorator."""
    return replica_training_mask


def replica_training_mask(exps_tr_masks, replica, experiments_index):
    """Save the boolean mask used to split data into training and validation
    for a given replica as a pandas DataFrame, indexed by
    :py:func:`validphys.results.experiments_index`. Can be used to reconstruct
    the training and validation data used in a fit.

    Parameters
    ----------
    exps_tr_masks: list[list[np.array]]
        Result of :py:func:`tr_masks` collected over experiments, which creates
        the nested structure. The outer list is
        len(group_dataset_inputs_by_experiment) and the inner-most list has an
        array for each dataset in that particular experiment - as defined by the
        metadata. The arrays should be 1-D boolean arrays which can be used as
        masks.
    replica: int
        The index of the replica.
    experiments_index: pd.MultiIndex
        Index returned by :py:func:`validphys.results.experiments_index`.


    Example
    -------
    >>> from validphys.api import API
    >>> ds_inp = [
    ...     {'dataset': 'NMC', 'frac': 0.75},
    ...     {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD'], 'frac': 0.75},
    ...     {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10, 'frac': 0.75}
    ... ]
    >>> API.replica_training_mask(dataset_inputs=ds_inp, replica=1, trvlseed=123, theoryid=162, use_cuts="nocuts", mcseed=None, genrep=False)
                         replica 1
    group dataset    id
    NMC   NMC        0        True
                    1        True
                    2       False
                    3        True
                    4        True
    ...                        ...
    CMS   CMSZDIFF12 45       True
                    46       True
                    47       True
                    48      False
                    49       True

    [345 rows x 1 columns]
    """
    all_masks = np.concatenate([ds_mask for exp_masks in exps_tr_masks for ds_mask in exp_masks])
    return pd.DataFrame(all_masks, columns=[f"replica {replica}"], index=experiments_index)


replicas_training_mask = collect("replica_training_mask", ("replicas",))


@table
def training_mask_table(training_mask):
    """Same as ``training_mask`` but with a table decorator"""
    return training_mask


def training_mask(replicas_training_mask):
    """Save the boolean mask used to split data into training and validation
    for each replica as a pandas DataFrame, indexed by
    :py:func:`validphys.results.experiments_index`. Can be used to reconstruct
    the training and validation data used in a fit.

    Parameters
    ----------
    replicas_exps_tr_masks: list[list[list[np.array]]]
        Result of :py:func:`replica_tr_masks` collected over replicas

    Example
    -------
    >>> from validphys.api import API
    >>> from reportengine.namespaces import NSList
    >>> # create namespace list for collects over replicas.
    >>> reps = NSList(list(range(1, 4)), nskey="replica")
    >>> ds_inp = [
    ...     {'dataset': 'NMC', 'frac': 0.75},
    ...     {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD'], 'frac': 0.75},
    ...     {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10, 'frac': 0.75}
    ... ]
    >>> API.training_mask(dataset_inputs=ds_inp, replicas=reps, trvlseed=123, theoryid=162, use_cuts="nocuts", mcseed=None, genrep=False)
                        replica 1  replica 2  replica 3
    group dataset    id
    NMC   NMC        0        True      False      False
                    1        True       True       True
                    2       False       True       True
                    3        True       True      False
                    4        True       True       True
    ...                        ...        ...        ...
    CMS   CMSZDIFF12 45       True       True       True
                    46       True      False       True
                    47       True       True       True
                    48      False       True       True
                    49       True       True       True

    [345 rows x 3 columns]

    """
    return pd.concat(replicas_training_mask, axis=1)


def _fitting_lagrange_dict(lambdadataset):
    """Loads a generic lambda dataset, often used for positivity and integrability datasets
    For more information see :py:func:`validphys.n3fit_data_utils.positivity_reader`.

    Parameters
    ----------
    lambdadataset: validphys.core.LagrangeSetSpec
        Positivity (or integrability) set which is to be loaded.

    Examples
    --------
    >>> from validphys.api import API
    >>> posdataset = {"dataset": "POSF2U", "maxlambda": 1e6}
    >>> pos = API.fitting_pos_dict(posdataset=posdataset, theoryid=162)
    >>> len(pos)
    9
    """
    integrability = isinstance(lambdadataset, IntegrabilitySetSpec)
    mode = "integrability" if integrability else "positivity"
    log.info("Loading %s dataset %s", mode, lambdadataset)
    positivity_datasets = validphys_group_extractor([lambdadataset], [])
    ndata = positivity_datasets[0].ndata
    return {
        "datasets": positivity_datasets,
        "trmask": np.ones(ndata, dtype=np.bool),
        "name": lambdadataset.name,
        "expdata": np.zeros((1, ndata)),
        "ndata": ndata,
        "positivity": True,
        "lambda": lambdadataset.maxlambda,
        "count_chi2": False,
        "integrability": integrability,
    }


def posdatasets_fitting_pos_dict(posdatasets=None):
    """Loads all positivity datasets. It is not allowed to be empty.

    Parameters
    ----------
    integdatasets: list[validphys.core.PositivitySetSpec]
        list containing the settings for the positivity sets. Examples of
        these can be found in the runcards located in n3fit/runcards. They have
        a format similar to ``dataset_input``.
    """
    if posdatasets is not None:
        return [_fitting_lagrange_dict(i) for i in posdatasets]
    log.warning("Not using any positivity datasets.")
    return None


# can't use collect here because integdatasets might not exist.
def integdatasets_fitting_integ_dict(integdatasets=None):
    """Loads the integrability datasets. Calls same function as
    :py:func:`fitting_pos_dict`, except on each element of
    ``integdatasets`` if ``integdatasets`` is not None.

    Parameters
    ----------
    integdatasets: list[validphys.core.IntegrabilitySetSpec]
        list containing the settings for the integrability sets. Examples of
        these can be found in the runcards located in n3fit/runcards. They have
        a format similar to ``dataset_input``.

    Examples
    --------
    >>> from validphys.api import API
    >>> integdatasets = [{"dataset": "INTEGXT3", "maxlambda": 1e2}]
    >>> res = API.integdatasets_fitting_integ_dict(integdatasets=integdatasets, theoryid=53)
    >>> len(res), len(res[0])
    (1, 9)
    >>> res = API.integdatasets_fitting_integ_dict(integdatasets=None)
    >>> print(res)
    None

    """
    if integdatasets is not None:
        return [_fitting_lagrange_dict(i) for i in integdatasets]
    log.warning("Not using any integrability datasets.")
    return None
