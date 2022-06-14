# -*- coding: utf-8 -*-
"""
Tools to obtain and analyse the pseudodata that was seen by the neural
networks during the fitting.
"""
from collections import namedtuple
import logging
import hashlib

import numpy as np
import pandas as pd

from validphys.covmats import INTRA_DATASET_SYS_NAME

from reportengine import collect

FILE_PREFIX = "datacuts_theory_fitting_"

log = logging.getLogger(__name__)

DataTrValSpec = namedtuple('DataTrValSpec', ['pseudodata', 'tr_idx', 'val_idx'])

context_index = collect("groups_index", ("fitcontext",))
read_fit_pseudodata = collect('read_replica_pseudodata', ('fitreplicas', 'fitcontextwithcuts'))
read_pdf_pseudodata = collect('read_replica_pseudodata', ('pdfreplicas', 'fitcontextwithcuts'))

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
    sorted_index = context_index.sortlevel(level=range(1,3))[0]

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
        val = pd.read_csv(validation_path, index_col=[0, 1, 2], sep="\t", names=[f"replica {replica}"])
    except FileNotFoundError as e:
        raise FileNotFoundError(
            "Could not find saved training and validation data files. "
            f"Please ensure {fit} was generated with the savepseudodata flag set to true"
        ) from e
    tr["type"], val["type"] = "training", "validation"

    pseudodata = pd.concat((tr, val))
    pseudodata.sort_index(level=range(1,3), inplace=True)

    pseudodata.index = sorted_index

    tr = pseudodata[pseudodata["type"]=="training"]
    val = pseudodata[pseudodata["type"]=="validation"]

    return DataTrValSpec(pseudodata.drop("type", axis=1), tr.index, val.index)


def make_replica(groups_dataset_inputs_loaded_cd_with_cuts, replica_mcseed, genrep=True):
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

    seed: int, None
        Seed used to initialise the numpy random number generator. If ``None`` then a random seed is
        allocated using the default numpy behaviour.

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
        return np.concatenate([cd.central_values for cd in groups_dataset_inputs_loaded_cd_with_cuts])

    # Seed the numpy RNG with the seed and the name of the datasets in this run
    name_salt = "-".join(i.setname for i in groups_dataset_inputs_loaded_cd_with_cuts)
    name_seed = int(hashlib.sha256(name_salt.encode()).hexdigest(), 16) % 10 ** 8
    rng = np.random.default_rng(seed=replica_mcseed+name_seed)

    # The inner while True loop is for ensuring a positive definite
    # pseudodata replica
    while True:
        pseudodatas = []
        special_add = []
        special_mult = []
        mult_shifts = []
        check_positive_masks = []
        for cd in groups_dataset_inputs_loaded_cd_with_cuts:
            # copy here to avoid mutating the central values.
            pseudodata = cd.central_values.to_numpy(copy=True)

            # add contribution from statistical uncertainty
            pseudodata += (cd.stat_errors.to_numpy() * rng.normal(size=cd.ndata))

            # ~~~ ADDITIVE ERRORS  ~~~
            add_errors = cd.additive_errors
            add_uncorr_errors = add_errors.loc[:, add_errors.columns=="UNCORR"].to_numpy()

            pseudodata += (add_uncorr_errors * rng.normal(size=add_uncorr_errors.shape)).sum(axis=1)

            # correlated within dataset
            add_corr_errors = add_errors.loc[:, add_errors.columns == "CORR"].to_numpy()
            pseudodata += add_corr_errors @ rng.normal(size=add_corr_errors.shape[1])

            # append the partially shifted pseudodata
            pseudodatas.append(pseudodata)
            # store the additive errors with correlations between datasets for later use
            special_add.append(
                add_errors.loc[:, ~add_errors.columns.isin(INTRA_DATASET_SYS_NAME)]
            )
            # ~~~ MULTIPLICATIVE ERRORS ~~~
            mult_errors = cd.multiplicative_errors
            mult_uncorr_errors = mult_errors.loc[:, mult_errors.columns == "UNCORR"].to_numpy()
            # convert to from percent to fraction
            mult_shift = (
                1 + mult_uncorr_errors * rng.normal(size=mult_uncorr_errors.shape) / 100
            ).prod(axis=1)

            mult_corr_errors = mult_errors.loc[:, mult_errors.columns == "CORR"].to_numpy()
            mult_shift *= (
                1 + mult_corr_errors * rng.normal(size=(1, mult_corr_errors.shape[1])) / 100
            ).prod(axis=1)

            mult_shifts.append(mult_shift)

            # store the multiplicative errors with correlations between datasets for later use
            special_mult.append(
                mult_errors.loc[:, ~mult_errors.columns.isin(INTRA_DATASET_SYS_NAME)]
            )

            # mask out the data we want to check are all positive
            if "ASY" in cd.commondataproc:
                check_positive_masks.append(np.zeros_like(pseudodata, dtype=bool))
            else:
                check_positive_masks.append(np.ones_like(pseudodata, dtype=bool))

        # non-overlapping systematics are set to NaN by concat, fill with 0 instead.
        special_add_errors = pd.concat(special_add, axis=0, sort=True).fillna(0).to_numpy()
        special_mult_errors = pd.concat(special_mult, axis=0, sort=True).fillna(0).to_numpy()


        all_pseudodata = (
            np.concatenate(pseudodatas, axis=0)
            + special_add_errors @ rng.normal(size=special_add_errors.shape[1])
        ) * (
            np.concatenate(mult_shifts, axis=0)
            * (1 + special_mult_errors * rng.normal(size=(1, special_mult_errors.shape[1])) / 100).prod(axis=1)
        )

        if np.all(all_pseudodata[np.concatenate(check_positive_masks, axis=0)] >= 0):
            break

    return all_pseudodata


def indexed_make_replica(groups_index, make_replica):
    """Index the make_replica pseudodata appropriately
    """

    return pd.DataFrame(make_replica, index=groups_index, columns=["data"])


_group_recreate_pseudodata = collect('indexed_make_replica', ('group_dataset_inputs_by_experiment',))
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