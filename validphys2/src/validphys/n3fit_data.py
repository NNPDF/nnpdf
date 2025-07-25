"""
n3fit_data.py

Providers which prepare the data ready for
:py:func:`n3fit.performfit.performfit`.
"""

from collections import abc, defaultdict
from copy import copy
import functools
import hashlib
import logging

import numpy as np
import pandas as pd

from reportengine import collect, namespaces
from reportengine.table import table
from validphys.core import IntegrabilitySetSpec, TupleComp
from validphys.n3fit_data_utils import validphys_group_extractor

log = logging.getLogger(__name__)


class Hashrray(TupleComp):
    """Wrapper class to hash a numpy array so it can be cached."""

    def __init__(self, array):
        self.array = array
        super().__init__(hash(self.array.tobytes()))


def _per_replica(f):
    """Decorator to be used on top of reportengine's decorators.
    It replaces the preparation step of the decorator with a custom function,
    which modifies the output behaviour when there is a collection of replicas.

    If there is no ``replica_path`` in the environment or collection over replicas
    this function does nothing. Otherwise, it removes the replica number from the
    output file and directs the output to ``replica_<replica>`` instead.
    """
    original_prepare = f.prepare

    def prepare_replica_path(*, spec, namespace, environment, **kwargs):
        if not hasattr(environment, "replica_path") or "replicas" not in namespace:
            return original_prepare(spec=spec, namespace=namespace, environment=environment)

        if not isinstance(namespace["replicas"], abc.Collection):
            return original_prepare(spec=spec, namespace=namespace, environment=environment)

        # Loop over the function input arguments to find the collection of replicas
        # pass down all other arguments unchanged
        rnumber = None
        new_nsspec = []
        for farg in spec.nsspec:
            if isinstance(farg, abc.Collection) and farg[0] == "replicas":
                rnumber = namespaces.value_from_spcec_ele(namespace, farg)
            else:
                new_nsspec.append(farg)
        if rnumber is None:
            raise ValueError("Wrong call to @_replica_table, no replica number found.")

        replica_path = environment.replica_path / f"replica_{rnumber}"

        new_env = copy(environment)
        new_env.table_folder = replica_path
        new_spec = spec._replace(nsspec=tuple(new_nsspec))

        return original_prepare(spec=new_spec, namespace=namespace, environment=new_env)

    f.prepare = prepare_replica_path

    return f


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


def replica_luxseed(replica, luxseed):
    """Generate the ``luxseed`` for a ``replica``.
    Identical to replica_nnseed but used for a different purpose.
    """
    return replica_nnseed(replica, luxseed)


class _Masks(TupleComp):
    """Class holding the training validation mask for a group of datasets
    If the same group of dataset receives the same trvlseed then the mask
    will be the same.
    This class holds said information so it can be reused easily, i.e.,
    ``group_name`` and ``seed`` define the ``masks``.

    In case of the diagonal basis, the ``rotation`` is a tuple of the eigenvalues and
    eigenvectors of the fitting covariance matrix.
    """

    def __init__(
        self,
        group_name,
        seed,
        tr_masks,
        vl_masks,
        diagonal_basis=False,
        eig_vals=None,
        diagonal_rotation=None,
    ):
        """
        Initialize the _Masks object.

        Parameters
        ----------
        group_name : str
            The name of the group of datasets.
        seed : int
            The seed used for generating the masks.
        tr_masks : list[np.array]
            List of boolean arrays representing the training masks.
        vl_masks : list[np.array]
            List of boolean arrays representing the validation masks.
        diagonal_basis : bool, optional
            Whether the masks are in the diagonal basis. Default is False.
        eig_vals : np.array, optional
            Eigenvalues of the covariance matrix, required if diagonal_basis is True.
        diagonal_rotation : np.array, optional
            Eigenvectors of the correlation matrix, required if diagonal_basis is True.
        """

        self.tr_masks = tr_masks
        self.vl_masks = vl_masks
        if diagonal_basis:
            self.eig_vals = eig_vals
            self.diagonal_rotation = diagonal_rotation

        super().__init__(group_name, seed)


def diagonal_masks(
    data, replica_trvlseed, dataset_inputs_fitting_covmat, diagonal_frac=1.0, threshold_eigvals=0
):

    # diagonalise the covariance matrix, eigenvalues appear in ascending order
    covmat = dataset_inputs_fitting_covmat

    # convert covmat to correlation
    diag_inv_sqrt = 1 / np.sqrt(np.diag(covmat))
    cormat = np.einsum("i, ij, j -> ij", diag_inv_sqrt, covmat, diag_inv_sqrt)

    # diagonalise the correlation matrix
    eig_vals, u_trans = np.linalg.eigh(cormat)
    u_trans = np.einsum("i, ik -> ik", diag_inv_sqrt, u_trans)
    ndata = len(eig_vals)

    # construct training mask by selecting a fraction of the eigenvalues
    tr_mask = np.random.random(ndata) < diagonal_frac
    vl_mask = ~tr_mask

    # discard the eigenvalues below the set threshold
    tr_mask[eig_vals < threshold_eigvals] = False
    vl_mask[eig_vals < threshold_eigvals] = False
    return _Masks(
        str(data),
        replica_trvlseed,
        [tr_mask],
        [vl_mask],
        diagonal_basis=True,
        eig_vals=eig_vals,
        diagonal_rotation=u_trans.T,
    )


def standard_masks(data, replica_trvlseed):
    """Generate the boolean masks used to split data into training and
    validation points. Returns a list of 1-D boolean arrays, one for each
    dataset. Each array has length equal to N_data, the datapoints which
    will be included in the training are ``True`` such that

        tr_data = data[tr_mask]

    """
    nameseed = int(hashlib.sha256(str(data).encode()).hexdigest(), 16) % 10**8
    nameseed += replica_trvlseed
    # TODO: update this to new random infrastructure.
    rng = np.random.Generator(np.random.PCG64(nameseed))

    trmask_partial = []
    vlmask_partial = []
    nomasking = True
    for dataset in data.datasets:
        # TODO: python commondata will not require this rubbish.
        # all data if cuts are None
        cuts = dataset.cuts
        ndata = len(cuts.load()) if cuts else dataset.commondata.ndata
        if ndata == 0:
            tr_mask = []
            vl_mask = []
            continue

        frac = dataset.frac
        # nomasking turns to False as soon as one frac is not equal to 1
        nomasking &= frac == 1.0
        # We do this so that a given dataset will always have the same number of points masked
        trmax = int(ndata * frac)
        if trmax == 0:
            # If that number is 0, then get 1 point with probability frac
            trmax = int(rng.random() < frac)
        tr_mask = np.concatenate([np.ones(trmax, dtype=bool), np.zeros(ndata - trmax, dtype=bool)])
        rng.shuffle(tr_mask)
        vl_mask = ~tr_mask
        trmask_partial.append(tr_mask)
        vlmask_partial.append(vl_mask)
    # if we are not masking, remove the seed from the object
    if nomasking:
        replica_trvlseed = None
    return _Masks(str(data), replica_trvlseed, trmask_partial, vlmask_partial)


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
                    mask.append(np.zeros(ndata, dtype=bool))
                # otherwise of ones
                else:
                    mask.append(np.ones(ndata, dtype=bool))
            list_folds.append(np.concatenate(mask))
    return list_folds


@functools.lru_cache
def fittable_datasets_masked(data):
    """Generate a list of :py:class:`validphys.n3fit_data_utils.FittableDataSet`
    from a group of dataset and the corresponding training/validation masks
    """
    # This is separated from fitting_data_dict so that we can cache the result
    # when the trvlseed is the same for all replicas (great for parallel replicas)
    return validphys_group_extractor(data.datasets)


def _hashed_dataset_inputs_fitting_covmat(dataset_inputs_fitting_covmat) -> Hashrray:
    """Wrap the covmat into a Hashrray for caches to work"""
    return Hashrray(dataset_inputs_fitting_covmat)


@functools.lru_cache
def _inv_covmat_prepared(masks, _hashed_dataset_inputs_fitting_covmat, diagonal_basis=False):
    """Returns the inverse covmats for training, validation and total
    attending to the right masks and whether it is diagonal or not.

    Since the masks and number of datapoints need to be treated for 1-point datasets
    it also returns the right ndata and masks for training and validation:

    inv_total, inv_training, inv_validation, ndata_tr, ndata_vl, mask_tr, mask_vl, diagonal_rotation
    """
    covmat = _hashed_dataset_inputs_fitting_covmat.array
    inv_total = np.linalg.inv(covmat)
    diagonal_rotation = None

    if diagonal_basis:
        log.info("working in diagonal basis.")

        # get the eigenvalues of the fit cormat (in ascending order)
        eig_vals = masks.eig_vals

        # rotate the experimental data to the diagonal basis of the cormat and obtain training/validation masks
        diagonal_rotation = masks.diagonal_rotation
        tr_mask = masks.tr_masks[0]
        vl_mask = masks.vl_masks[0]

        # apply the training/validation masks to the eigenvalues and take the inverse
        # this does not give the inverse of the covmat as the variable name might suggest,
        # but we call it this way anyway as this needs to be returned at the end
        invcovmat_tr = np.diag(1 / eig_vals[tr_mask])
        invcovmat_vl = np.diag(1 / eig_vals[vl_mask])

        # obtain the number of data points in the training/validation sets
        ndata_tr = invcovmat_tr.shape[0]
        ndata_vl = invcovmat_vl.shape[0]

    else:
        # In the fittable datasets the fktables masked for 1-point datasets will be set to 0
        # Here we want to have the data both in training and validation,
        # but set to 0 the data, so that it doesn't affect the chi2 value.

        zero_tr = []
        zero_vl = []
        idx = 0
        for data_mask in masks.tr_masks:
            dlen = len(data_mask)
            if dlen == 1:
                if data_mask[0]:
                    zero_vl.append(idx)
                else:
                    zero_tr.append(idx)
            idx += dlen

        tr_mask = np.concatenate(masks.tr_masks)
        vl_mask = ~tr_mask

        # Now set to true the masks
        tr_mask[zero_tr] = True
        vl_mask[zero_vl] = True
        # And prepare the index to 0 the (inverse) covmat
        data_zero_tr = np.cumsum(tr_mask)[zero_tr] - 1
        data_zero_vl = np.cumsum(vl_mask)[zero_vl] - 1

        covmat_tr = covmat[tr_mask].T[tr_mask]
        covmat_vl = covmat[vl_mask].T[vl_mask]

        # Remove possible correlations for 1-point datasets that should've been masked out
        covmat_tr[data_zero_tr, :] = covmat_tr[:, data_zero_tr] = 0.0
        covmat_vl[data_zero_vl, :] = covmat_vl[:, data_zero_vl] = 0.0
        # Set the diagonal to 1 to avoid infinities or inconsistencies when computing the inverse
        covmat_tr[data_zero_tr, data_zero_tr] = 1.0
        covmat_vl[data_zero_vl, data_zero_vl] = 1.0

        diag_inv_sqrt_covmat_tr = 1 / np.sqrt(np.diag(covmat_tr))
        diag_inv_sqrt_covmat_vl = 1 / np.sqrt(np.diag(covmat_vl))
        cormat_tr = np.einsum(
            "i, ij, j -> ij", diag_inv_sqrt_covmat_tr, covmat_tr, diag_inv_sqrt_covmat_tr
        )
        cormat_vl = np.einsum(
            "i, ij, j -> ij", diag_inv_sqrt_covmat_vl, covmat_vl, diag_inv_sqrt_covmat_vl
        )
        invcovmat_tr = np.einsum(
            "i, ij, j -> ij",
            diag_inv_sqrt_covmat_tr,
            np.linalg.inv(cormat_tr),
            diag_inv_sqrt_covmat_tr,
        )
        invcovmat_vl = np.einsum(
            "i, ij, j -> ij",
            diag_inv_sqrt_covmat_vl,
            np.linalg.inv(cormat_vl),
            diag_inv_sqrt_covmat_vl,
        )

        # Set to 0 the points in the diagonal that were left as 1
        invcovmat_tr[np.ix_(data_zero_tr, data_zero_tr)] = 0.0
        invcovmat_vl[np.ix_(data_zero_vl, data_zero_vl)] = 0.0

        ndata_tr = np.count_nonzero(tr_mask)
        ndata_vl = np.count_nonzero(vl_mask)

        # And subtract them for ndata
        ndata_tr -= len(data_zero_tr)
        ndata_vl -= len(data_zero_vl)

    return (
        inv_total,
        invcovmat_tr,
        invcovmat_vl,
        ndata_tr,
        ndata_vl,
        tr_mask,
        vl_mask,
        diagonal_rotation,
    )


def fitting_data_dict(
    data,
    make_replica,
    dataset_inputs_loaded_cd_with_cuts,
    dataset_inputs_fitting_covmat,
    _inv_covmat_prepared,
    kfold_masks,
    fittable_datasets_masked,
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
    fittable_datasets = fittable_datasets_masked

    inv_true, invcovmat_tr, invcovmat_vl, ndata_tr, ndata_vl, tr_mask, vl_mask, diag_rot = (
        _inv_covmat_prepared
    )

    if diag_rot is not None:
        expdata = diag_rot @ expdata

    expdata_tr = expdata[tr_mask].reshape(1, -1)
    expdata_vl = expdata[vl_mask].reshape(1, -1)

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
        "covmat": dataset_inputs_fitting_covmat,
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
        "data_transformation": diag_rot,
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


replicas_training_pseudodata = collect("training_pseudodata", ("replicas",))
replicas_validation_pseudodata = collect("validation_pseudodata", ("replicas",))
replicas_pseudodata = collect("pseudodata_table", ("replicas",))
replicas_nnseed_fitting_data_dict = collect("replica_nnseed_fitting_data_dict", ("replicas",))
groups_replicas_indexed_make_replica = collect(
    "indexed_make_replica", ("replicas", "group_dataset_inputs_by_metadata")
)

experiment_indexed_make_replica = collect(
    "indexed_make_replica", ("group_dataset_inputs_by_metadata",)
)


def replica_pseudodata(experiment_indexed_make_replica, replica):
    """Creates a pandas DataFrame containing the generated pseudodata.
    The index is :py:func:`validphys.results.experiments_index` and the columns
    is the replica numbers.

    Notes
    -----
    Whilst running ``n3fit``, this action will only be called if
    `fitting::savepseudodata` is `true` (as per the default setting)
    The table can be found in the replica folder i.e. <fit dir>/nnfit/replica_*/
    """
    df = pd.concat(experiment_indexed_make_replica)
    df.columns = [f"replica {replica}"]
    return df


@_per_replica
@table
def pseudodata_table(replica_pseudodata):
    """Save the pseudodata for the given replica.
    Deactivate by setting ``fitting::savepseudodata: False``
    from within the fit runcard.
    """
    return replica_pseudodata


@_per_replica
@table
def training_pseudodata(replica_pseudodata, replica_mask):
    """Save the training data for the given replica.
    Deactivate by setting ``fitting::savepseudodata: False``
    from within the fit runcard.

    See Also
    --------
    :py:func:`validphys.n3fit_data.validation_pseudodata`
    """
    return replica_pseudodata.loc[replica_mask[0].values]


@_per_replica
@table
def validation_pseudodata(replica_pseudodata, replica_mask):
    """Save the training data for the given replica.
    Deactivate by setting ``fitting::savepseudodata: False``
    from within the fit runcard.

    See Also
    --------
    :py:func:`validphys.n3fit_data.validation_pseudodata`
    """
    return replica_pseudodata.loc[replica_mask[1].values]


exps_masks = collect("masks", ("group_dataset_inputs_by_metadata",))
replicas_exps_masks = collect("exps_masks", ("replicas",))


@table
def replica_mask_table(replica_mask):
    """Same as ``replica_training_mask`` but with a table decorator."""
    return replica_mask


def replica_mask(exps_masks, replica, experiments_index, diagonal_basis=False):
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
    ...     {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy', 'frac': 0.75},
    ...     {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy_theory', 'frac': 0.75},
    ...     {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac':('NRM',), 'frac': 0.75},
    ... ]
    >>> API.replica_training_mask(dataset_inputs=ds_inp, replica=1, trvlseed=123, theoryid=40_000_000, use_cuts="nocuts", mcseed=None, genrep=False)
                                        replica 1
    group dataset                       id
    NMC   NMC_NC_NOTFIXED_P_EM-SIGMARED 0        True
                                        1        True
                                        2        True
                                        3        True
                                        4       False
    ...                                           ...
    CMS   CMS_Z0J_8TEV_PT-Y             45       True
                                        46       True
                                        47       True
                                        48       True
                                        49       True

    [343 rows x 1 columns]
    """

    all_tr_masks = np.concatenate(
        [ds_mask for exp_masks in exps_masks for ds_mask in exp_masks.tr_masks]
    )
    all_vl_masks = np.concatenate(
        [ds_mask for exp_masks in exps_masks for ds_mask in exp_masks.vl_masks]
    )

    index = (
        [f"eigenmode {i}" for i in range(len(all_tr_masks))]
        if diagonal_basis
        else experiments_index
    )

    df_tr = pd.DataFrame(all_tr_masks, columns=[f"replica {replica}"], index=index)
    df_vl = pd.DataFrame(all_vl_masks, columns=[f"replica {replica}"], index=index)

    return df_tr, df_vl


def replica_validation_mask(exps_tr_masks, replica, experiments_index, diagonal_basis=False):
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

    all_masks = np.concatenate(
        [ds_mask for exp_masks in exps_masks for ds_mask.vl_masks in exp_masks]
    )
    if diagonal_basis:
        return pd.DataFrame(
            all_masks,
            columns=[f"replica {replica}"],
            index=[f"eigenmode {i}" for i in range(len(all_masks))],
        )
    else:
        return pd.DataFrame(all_masks, columns=[f"replica {replica}"], index=experiments_index)


replicas_mask = collect("replica_mask", ("replicas",))


@table
def training_mask_table(training_mask):
    """Same as ``training_mask`` but with a table decorator"""
    return training_mask


def training_mask(replicas_mask):
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
    ...     {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy', 'frac': 0.75},
    ...     {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy_theory', 'frac': 0.75},
    ...     {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac':('NRM',), 'frac': 0.75},
    ... ]
    >>> API.training_mask(dataset_inputs=ds_inp, nreplica=3, trvlseed=123, theoryid=40_000_000, use_cuts="nocuts", mcseed=None, genrep=False)
                                                    replica 1  replica 2  replica 3
        group dataset                       id
        NMC   NMC_NC_NOTFIXED_P_EM-SIGMARED 0        True      False      False
                                            1        True       True       True
                                            2        True      False       True
                                            3        True       True      False
                                            4       False       True       True
        ...                                           ...        ...        ...
        CMS   CMS_Z0J_8TEV_PT-Y             45       True      False       True
                                            46       True       True       True
                                            47       True      False       True
                                            48       True       True       True
                                            49       True      False       True
        [343 rows x 3 columns]

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
    >>> posdataset = {"dataset": "NNPDF_POS_2P24GEV_F2U", "maxlambda": 1e6}
    >>> pos = API.fitting_pos_dict(posdataset=posdataset, theoryid=40_000_000)
    >>> len(pos)
    9
    """
    integrability = isinstance(lambdadataset, IntegrabilitySetSpec)
    mode = "integrability" if integrability else "positivity"
    log.info("Loading %s dataset %s", mode, lambdadataset)
    positivity_datasets = validphys_group_extractor([lambdadataset])
    ndata = positivity_datasets[0].ndata
    return {
        "datasets": positivity_datasets,
        "trmask": np.ones(ndata, dtype=bool),
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
