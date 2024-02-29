# -*- coding: utf-8 -*-
"""
Tools to obtain and analyse the pseudodata that was seen by the neural
networks during the fitting.
"""
from collections import namedtuple
import hashlib
import logging

import numpy as np
import pandas as pd

from reportengine import collect
from validphys.covmats import (
    INTRA_DATASET_SYS_NAME,
    dataset_inputs_covmat_from_systematics,
    sqrt_covmat,
)

FILE_PREFIX = "datacuts_theory_fitting_"

log = logging.getLogger(__name__)

DataTrValSpec = namedtuple('DataTrValSpec', ['pseudodata', 'tr_idx', 'val_idx'])

context_index = collect("groups_index", ("fitcontext",))
read_fit_pseudodata = collect('read_replica_pseudodata', ('fitreplicas', 'fitcontextwithcuts'))
read_pdf_pseudodata = collect('read_replica_pseudodata', ('pdfreplicas', 'fitcontextwithcuts'))


class ReplicaGenerationError(Exception):
    pass


def read_replica_pseudodata(fit, context_index, replica):
    """Function to handle the reading of training and validation splits for a fit that has been
    produced with the ``savepseudodata`` flag set to ``True``.

    The data is read from the PDF to handle the mixing introduced by ``postfit``.

    The data files are concatenated to return all the data that went into a fit. The training and validation
    indices are also returned so one can access the splits using pandas indexing.

    Raises
    ------
    FileNotFoundError
        If the training or validation files for the PDF set cannot be found.
    CheckError
        If the ``use_cuts`` flag is not set to ``fromfit``

    Returns
    -------
    data_indices_list: list[namedtuple]
        List of ``namedtuple`` where each entry corresponds to a given replica. Each element contains
        attributes ``pseudodata``, ``tr_idx``, and ``val_idx``. The latter two being used to slice
        the former to return training and validation data respectively.

    Example
    -------
    >>> from validphys.api import API
    >>> data_indices_list = API.read_fit_pseudodata(fit="pseudodata_test_fit_n3fit")
    >>> len(data_indices_list) # Same as nrep
    10
    >>> rep_info = data_indices_list[0]
    >>> rep_info.pseudodata.loc[rep_info.tr_idx].head()
                                replica 1
    group dataset           id
    ATLAS ATLASZPT8TEVMDIST 1   30.665835
                            3   15.795880
                            4    8.769734
                            5    3.117819
                            6    0.771079
    """
    # List of length 1 due to the collect
    context_index = context_index[0]
    # The [0] is because of how pandas handles sorting a MultiIndex
    sorted_index = context_index.sortlevel(level=range(1, 3))[0]

    log.debug(f"Reading pseudodata & training/validation splits from {fit.name}.")
    replica_path = fit.path / "nnfit" / f"replica_{replica}"

    training_path = replica_path / (FILE_PREFIX + "training_pseudodata.csv")
    validation_path = replica_path / (FILE_PREFIX + "validation_pseudodata.csv")

    try:
        tr = pd.read_csv(training_path, index_col=[0, 1, 2], sep="\t", header=0)
        val = pd.read_csv(validation_path, index_col=[0, 1, 2], sep="\t", header=0)
    except FileNotFoundError:
        # Old 3.1 style fits had pseudodata called training.dat and validation.dat
        training_path = replica_path / "training.dat"
        validation_path = replica_path / "validation.dat"
        tr = pd.read_csv(training_path, index_col=[0, 1, 2], sep="\t", names=[f"replica {replica}"])
        val = pd.read_csv(
            validation_path, index_col=[0, 1, 2], sep="\t", names=[f"replica {replica}"]
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(
            "Could not find saved training and validation data files. "
            f"Please ensure {fit} was generated with the savepseudodata flag set to true"
        ) from e
    tr["type"], val["type"] = "training", "validation"

    pseudodata = pd.concat((tr, val))
    pseudodata.sort_index(level=range(1, 3), inplace=True)

    pseudodata.index = sorted_index

    tr = pseudodata[pseudodata["type"] == "training"]
    val = pseudodata[pseudodata["type"] == "validation"]

    return DataTrValSpec(pseudodata.drop("type", axis=1), tr.index, val.index)


def make_replica(
    groups_dataset_inputs_loaded_cd_with_cuts,
    replica_mcseed,
    dataset_inputs_sampling_covmat,
    sep_mult,
    genrep=True,
    max_tries=int(1e6),
):
    """Function that takes in a list of :py:class:`validphys.coredata.CommonData`
    objects and returns a pseudodata replica accounting for
    possible correlations between systematic uncertainties.

    The function loops until positive definite pseudodata is generated for any
    non-asymmetry datasets. In the case of an asymmetry dataset negative values are
    permitted so the loop block executes only once.

    Parameters
    ---------
    groups_dataset_inputs_loaded_cd_with_cuts: list[:py:class:`validphys.coredata.CommonData`]
        List of CommonData objects which stores information about systematic errors,
        their treatment and description, for each dataset.

    replica_mcseed: int, None
        Seed used to initialise the numpy random number generator. If ``None`` then a random seed is
        allocated using the default numpy behaviour.

    dataset_inputs_sampling_covmat: np.array
        Full covmat to be used. It can be either only experimental or also theoretical.

    separate_multiplicative: bool
        Specifies whether computing the shifts with the full covmat or separating multiplicative
        errors (in the latter case remember to generate the covmat coherently)

    genrep: bool
        Specifies whether computing replicas or not

    max_tries: int
        The stochastic nature of replica generation means one can obtain (unphysical) negative predictions.
        If after max_tries (default=1e6) no physical configuration is found,
        it will raise a :py:class:`ReplicaGenerationError`

    Returns
    -------
    pseudodata: np.array
        Numpy array which is N_dat (where N_dat is the combined number of data points after cuts)
        containing monte carlo samples of data centered around the data central value.

    Example
    -------
    >>> from validphys.api import API
    >>> pseudodata = API.make_replica(
                                    dataset_inputs=[{"dataset":"NMC"}, {"dataset": "NMCPD"}],
                                    use_cuts="nocuts",
                                    theoryid=53,
                                    replica=1,
                                    mcseed=123,
                                    genrep=True,
                                )
    array([0.25640033, 0.25986534, 0.27165461, 0.29001009, 0.30863588,
       0.30100351, 0.31781208, 0.30827054, 0.30258217, 0.32116842,
       0.34206012, 0.31866286, 0.2790856 , 0.33257621, 0.33680007,
    """
    if not genrep:
        return np.concatenate(
            [cd.central_values for cd in groups_dataset_inputs_loaded_cd_with_cuts]
        )

    # Seed the numpy RNG with the seed and the name of the datasets in this run
    name_salt = "-".join(i.setname for i in groups_dataset_inputs_loaded_cd_with_cuts)
    name_seed = int(hashlib.sha256(name_salt.encode()).hexdigest(), 16) % 10**8
    rng = np.random.default_rng(seed=replica_mcseed + name_seed)
    # construct covmat
    covmat = dataset_inputs_sampling_covmat
    covmat_sqrt = sqrt_covmat(covmat)
    # Loading the data
    pseudodatas = []
    check_positive_masks = []
    nonspecial_mult = []
    special_mult = []
    for cd in groups_dataset_inputs_loaded_cd_with_cuts:
        # copy here to avoid mutating the central values.
        pseudodata = cd.central_values.to_numpy()

        pseudodatas.append(pseudodata)
        # Separation of multiplicative errors. If separate_multiplicative is True also the exp_covmat is produced
        # without multiplicative errors
        if sep_mult:
            mult_errors = cd.multiplicative_errors
            mult_uncorr_errors = mult_errors.loc[:, mult_errors.columns == "UNCORR"].to_numpy()
            mult_corr_errors = mult_errors.loc[:, mult_errors.columns == "CORR"].to_numpy()
            nonspecial_mult.append((mult_uncorr_errors, mult_corr_errors))
            special_mult.append(
                mult_errors.loc[:, ~mult_errors.columns.isin(INTRA_DATASET_SYS_NAME)]
            )
        if "ASY" in cd.commondataproc:
            check_positive_masks.append(np.zeros_like(pseudodata, dtype=bool))
        else:
            check_positive_masks.append(np.ones_like(pseudodata, dtype=bool))
    # concatenating special multiplicative errors, pseudodatas and positive mask
    if sep_mult:
        special_mult_errors = pd.concat(special_mult, axis=0, sort=True).fillna(0).to_numpy()
    all_pseudodata = np.concatenate(pseudodatas, axis=0)
    full_mask = np.concatenate(check_positive_masks, axis=0)
    # The inner while True loop is for ensuring a positive definite
    # pseudodata replica
    for _ in range(max_tries):
        mult_shifts = []
        # Prepare the per-dataset multiplicative shifts
        for mult_uncorr_errors, mult_corr_errors in nonspecial_mult:
            # convert to from percent to fraction
            mult_shift = (
                1 + mult_uncorr_errors * rng.normal(size=mult_uncorr_errors.shape) / 100
            ).prod(axis=1)

            mult_shift *= (
                1 + mult_corr_errors * rng.normal(size=(1, mult_corr_errors.shape[1])) / 100
            ).prod(axis=1)

            mult_shifts.append(mult_shift)

        # If sep_mult is true then the multiplicative shifts were not included in the covmat
        shifts = covmat_sqrt @ rng.normal(size=covmat.shape[1])
        mult_part = 1.0
        if sep_mult:
            special_mult = (
                1 + special_mult_errors * rng.normal(size=(1, special_mult_errors.shape[1])) / 100
            ).prod(axis=1)
            mult_part = np.concatenate(mult_shifts, axis=0) * special_mult
        # Shifting pseudodata
        shifted_pseudodata = (all_pseudodata + shifts) * mult_part
        # positivity control
        if np.all(shifted_pseudodata[full_mask] >= 0):
            return shifted_pseudodata

    dfail = " ".join(i.setname for i in groups_dataset_inputs_loaded_cd_with_cuts)
    log.error(f"Error generating replicas for the group: {dfail}")
    raise ReplicaGenerationError(f"No valid replica found after {max_tries} attempts")


def indexed_make_replica(groups_index, make_replica):
    """Index the make_replica pseudodata appropriately"""

    return pd.DataFrame(make_replica, index=groups_index, columns=["data"])


def level0_commondata_wc(data, fakepdf):
    """
    Given a validphys.core.DataGroupSpec object, load commondata and
    generate a new commondata instance with central values replaced
    by fakepdf prediction

    Parameters
    ----------

    data : validphys.core.DataGroupSpec

    fakepdf: validphys.core.PDF

    Returns
    -------
    list
        list of validphys.coredata.CommonData instances corresponding to
        all datasets within one experiment. The central value is replaced
        by Level 0 fake data.

    Example
    -------
    >>> from validphys.api import API
    >>> API.level0_commondata_wc(dataset_inputs = [{"dataset":"NMC"}], use_cuts="internal", theoryid=200,fakepdf = "NNPDF40_nnlo_as_01180")

    [CommonData(setname='NMC', ndata=204, commondataproc='DIS_NCE', nkin=3, nsys=16)]
    """
    from validphys.covmats import dataset_t0_predictions

    level0_commondata_instances_wc = []

    # ==== Load validphys.coredata.CommonData instance with cuts ====#

    for dataset in data.datasets:
        commondata_wc = dataset.commondata.load_commondata_instance()
        if dataset.cuts is not None:
            cuts = dataset.cuts.load()
            commondata_wc = commondata_wc.with_cuts(cuts)

        # == Generate a new CommonData instance with central value given by Level 0 data generated with fakepdf ==#

        t0_prediction = dataset_t0_predictions(
            dataset=dataset, t0set=fakepdf
        )  # N.B. cuts already applied to th. pred.
        level0_commondata_instances_wc.append(commondata_wc.with_central_value(t0_prediction))

    return level0_commondata_instances_wc

def level0_commondata_wc_patched(data, fakepdf):
    import ipdb; ipdb.set_trace()
    return level0_commondata_wc(data, fakepdf)

def make_level1_data(data, level0_commondata_wc_patched, filterseed, data_index, sep_mult):
    """
    Given a list of Level 0 commondata instances, return the
    same list with central values replaced by Level 1 data.

    Level 1 data is generated using validphys.make_replica.
    The covariance matrix, from which the stochastic Level 1
    noise is sampled, is built from Level 0 commondata
    instances (level0_commondata_wc). This, in particular,
    means that the multiplicative systematics are generated
    from the Level 0 central values.

    Note that the covariance matrix used to generate Level 2
    pseudodata is consistent with the one used at Level 1
    up to corrections of the order eta * eps, where eta and
    eps are defined as shown below:

    Generate L1 data: L1 = L0 + eta, eta ~ N(0,CL0)
    Generate L2 data: L2_k = L1 + eps_k, eps_k ~ N(0,CL1)

    where CL0 and CL1 means that the multiplicative entries
    have been constructed from Level 0 and Level 1 central
    values respectively.


    Parameters
    ----------

    data : validphys.core.DataGroupSpec

    level0_commondata_wc : list
                        list of validphys.coredata.CommonData instances corresponding to
                        all datasets within one experiment. The central value is replaced
                        by Level 0 fake data. Cuts already applied.

    filterseed : int
                random seed used for the generation of Level 1 data

    data_index : pandas.MultiIndex

    Returns
    -------
    list
        list of validphys.coredata.CommonData instances corresponding to
        all datasets within one experiment. The central value is replaced
        by Level 1 fake data.

    Example
    -------

    >>> from validphys.api import API
    >>> dataset='NMC'
    >>> l1_cd = API.make_level1_data(dataset_inputs = [{"dataset":dataset}],use_cuts="internal", theoryid=200,
                             fakepdf = "NNPDF40_nnlo_as_01180",filterseed=1)
    >>> l1_cd
    [CommonData(setname='NMC', ndata=204, commondataproc='DIS_NCE', nkin=3, nsys=16)]
    """

    dataset_input_list = list(data.dsinputs)

    covmat = dataset_inputs_covmat_from_systematics(
        level0_commondata_wc_patched,
        dataset_input_list,
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=sep_mult,
    )

    # ================== generation of Level1 data ======================#
    level1_data = make_replica(
        level0_commondata_wc_patched, filterseed, covmat, sep_mult=sep_mult, genrep=True
    )

    indexed_level1_data = indexed_make_replica(data_index, level1_data)

    dataset_order = {cd.setname: i for i, cd in enumerate(level0_commondata_wc_patched)} 

    # ===== create commondata instances with central values given by pseudo_data =====#
    level1_commondata_dict = {c.setname: c for c in level0_commondata_wc_patched}
    level1_commondata_instances_wc = []

    for xx, grp in indexed_level1_data.groupby('dataset'):
        level1_commondata_instances_wc.append(
            level1_commondata_dict[xx].with_central_value(grp.values)
        )
    # sort back so as to mantain same order as in level0_commondata_wc
    level1_commondata_instances_wc.sort(key=lambda x: dataset_order[x.setname])
    
    return level1_commondata_instances_wc


_group_recreate_pseudodata = collect(
    'indexed_make_replica', ('group_dataset_inputs_by_experiment',)
)
_recreate_fit_pseudodata = collect('_group_recreate_pseudodata', ('fitreplicas', 'fitenvironment'))
_recreate_pdf_pseudodata = collect('_group_recreate_pseudodata', ('pdfreplicas', 'fitenvironment'))

fit_tr_masks = collect('replica_training_mask_table', ('fitreplicas', 'fitenvironment'))
pdf_tr_masks = collect('replica_training_mask_table', ('pdfreplicas', 'fitenvironment'))
make_replicas = collect('make_replica', ('replicas',))
fitted_make_replicas = collect('make_replica', ('pdfreplicas',))
indexed_make_replicas = collect('indexed_make_replica', ('replicas',))


def recreate_fit_pseudodata(_recreate_fit_pseudodata, fitreplicas, fit_tr_masks):
    """Function used to reconstruct the pseudodata seen by each of the
    Monte Carlo fit replicas.

    Returns
    -------
    res : list[namedtuple]
          List of namedtuples, each of which contains a dataframe
          containing all the data points, the training indices, and
          the validation indices.

    Example
    -------
    >>> from validphys.api import API
    >>> API.recreate_fit_pseudodata(fit="pseudodata_test_fit_n3fit")

    Notes
    -----
    - This function does not account for the postfit reshuffling.

    See Also
    --------
    :py:func:`validphys.pseudodata.recreate_pdf_pseudodata`
    """
    res = []
    for pseudodata, mask, rep in zip(_recreate_fit_pseudodata, fit_tr_masks, fitreplicas):
        df = pd.concat(pseudodata)
        df.columns = [f"replica {rep}"]
        tr_idx = df.loc[mask.values].index
        val_idx = df.loc[~mask.values].index
        res.append(DataTrValSpec(df, tr_idx, val_idx))
    return res


def recreate_pdf_pseudodata(_recreate_pdf_pseudodata, pdfreplicas, pdf_tr_masks):
    """Like :py:func:`validphys.pseudodata.recreate_fit_pseudodata`
    but accounts for the postfit reshuffling of replicas.

    Returns
    -------
    res : list[namedtuple]
          List of namedtuples, each of which contains a dataframe
          containing all the data points, the training indices, and
          the validation indices.

    Example
    -------
    >>> from validphys.api import API
    >>> API.recreate_pdf_pseudodata(fit="pseudodata_test_fit_n3fit")

    See Also
    --------
    :py:func:`validphys.pseudodata.recreate_fit_pseudodata`
    """
    return recreate_fit_pseudodata(_recreate_pdf_pseudodata, pdfreplicas, pdf_tr_masks)


pdf_tr_masks_no_table = collect('replica_training_mask', ('pdfreplicas', 'fitenvironment'))


def recreate_pdf_pseudodata_no_table(_recreate_pdf_pseudodata, pdfreplicas, pdf_tr_masks_no_table):
    return recreate_pdf_pseudodata(_recreate_pdf_pseudodata, pdfreplicas, pdf_tr_masks_no_table)
