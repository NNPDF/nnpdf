# -*- coding: utf-8 -*-
"""
Tools to obtain and analyse the pseudodata that was seen by the neural
networks during the fitting.
"""
from collections import namedtuple
import logging
import multiprocessing as mp
import os
import pathlib

import numpy as np
import pandas as pd

from validphys.checks import check_cuts_fromfit
from validphys.covmats import INTRA_DATASET_SYS_NAME

from reportengine import collect

import validphys.n3fit_data_utils as reader

log = logging.getLogger(__name__)

DataTrValSpec = namedtuple('DataTrValSpec', ['pseudodata', 'tr_idx', 'val_idx'])

context_index = collect("groups_index", ("fitcontext",))

@check_cuts_fromfit
def read_fit_pseudodata(fitcontext, context_index):
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
    >>> data_indices_list = API.read_fit_pseudodata(fit="NNPDF31_nnlo_as_0118_DISonly_pseudodata", use_cuts="fromfit")
    >>> len(data_indices_list) # Same as nrep
    100
    >>> rep_info = data_indices_list[0]
    >>> rep_info.pseudodata.loc[rep_info.tr_idx].head()
                          data
    group dataset id
    BCDMS BCDMSD  0   0.371510
                1   0.365659
                2   0.350234
                4   0.355560
                6   0.346234
                        data
    """
    # List of length 1 due to the collect
    context_index = context_index[0]
    # The [0] is because of how pandas handles sorting a MultiIndex
    sorted_index = context_index.sortlevel(level=range(1,3))[0]

    pdf = fitcontext["pdf"]
    log.debug(f"Using same pseudodata & training/validation splits as {pdf.name}.")
    nrep = len(pdf)
    path = pathlib.Path(pdf.infopath)

    data_indices_list = []
    for rep_number in range(1, nrep):
        # This is a symlink (usually).
        replica = path.with_name(pdf.name + "_" + str(rep_number).zfill(4) + ".dat")
        # we resolve the symlink
        if replica.parent.is_symlink():
            replica = pathlib.Path(os.path.realpath(replica))

        training_path = replica.with_name("training.dat")
        validation_path = replica.with_name("validation.dat")

        try:
            tr = pd.read_csv(training_path, index_col=[0, 1, 2], sep="\t", names=["data"])
            val = pd.read_csv(validation_path, index_col=[0, 1, 2], sep="\t", names=["data"])
        except FileNotFoundError as e:
            raise FileNotFoundError(
                "Could not find saved training and validation data files. "
                f"Please ensure {pdf} was generated with the savepseudodata flag set to true"
            ) from e
        tr["type"], val["type"] = "training", "validation"

        pseudodata = pd.concat((tr, val))
        pseudodata.sort_index(level=range(1,3), inplace=True)

        pseudodata.index = sorted_index

        tr = pseudodata[pseudodata["type"]=="training"]
        val = pseudodata[pseudodata["type"]=="validation"]

        data_indices_list.append(
            DataTrValSpec(pseudodata.drop("type", axis=1), tr.index, val.index)
        )

    return data_indices_list


def make_replica(dataset_inputs_loaded_cd_with_cuts, replica_mcseed):
    """Function that takes in a list of :py:class:`validphys.coredata.CommonData`
    objects and returns a pseudodata replica accounting for
    possible correlations between systematic uncertainties.

    The function loops until positive definite pseudodata is generated for any
    non-asymmetry datasets. In the case of an asymmetry dataset negative values are
    permitted so the loop block executes only once.

    Parameters
    ---------
    dataset_inputs_loaded_cd_with_cuts: list[:py:class:`validphys.coredata.CommonData`]
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
                                    theoryid=53
                                )
    array([0.25721162, 0.2709698 , 0.27525357, 0.28903442, 0.3114298 ,
        0.3005844 , 0.3184538 , 0.31094522, 0.30750703, 0.32673155,
        0.34843355, 0.34730928, 0.3090914 , 0.32825111, 0.3485292 ,
    """
    # Seed the numpy RNG with the seed.
    rng = np.random.default_rng(seed=replica_mcseed)

    # The inner while True loop is for ensuring a positive definite
    # pseudodata replica
    while True:
        pseudodatas = []
        special_add = []
        special_mult = []
        mult_shifts = []
        check_positive_masks = []
        for cd in dataset_inputs_loaded_cd_with_cuts:
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


def _make_replica_task(dataset_inputs, theoryid, fit, mcseed, genrep, API, replica):
    """A global function for
    :py:func:`validphys.pseudodata.indexed_make_replica`, to be used in
    conjunction with :py:func:`validphys.pseudodata.recreate_fit_pseudodata`.
    Defined at top level to be picklable for use with multiprocessing.
    """
    res = API.indexed_make_replica(
        dataset_inputs=dataset_inputs,
        theoryid=theoryid,
        use_cuts="fromfit",
        fit=str(fit),
        mcseed=mcseed,
        genrep=genrep,
        replica=replica,
    )
    return res


def recreate_fit_pseudodata(fit, fitinputcontext, num_fitted_replicas, NPROC=None):
    """Function that recreates the pseudodata seen by each Monte Carlo replica.
    Returns a dataframe with each column labelled ``replica i`` for ``i`` 1 through
    to the number of fitted replicas.

    # TODO: Finish this example
    Example
    -------

    Notes
    -----
    By default this function computes the replicas in parallel using the maximum
    number of available CPUs. Consider setting the ``NPROC`` flag to something
    smaller to leave resources available.
    """
    # The + 1 coming from the fact that we wish to
    # include the last replica
    replicas = range(1, num_fitted_replicas + 1)

    mcseed = fit.as_input()["mcseed"]
    genrep = fit.as_input()["genrep"]

    dataset_inputs = fit.as_input()['dataset_inputs']
    theoryid = fitinputcontext["theoryid"].id

    # Because of the fact that only global functions can be
    # pickled, we need to define the task as a high level function.
    # As such we need to do this annoying repeating of function arguments,
    # since only `replica` changes between monte carlo replicas, but
    # dataset_inputs, fit, mcseed and genrep stay identical. See:
    # https://stackoverflow.com/questions/5442910/how-to-use-multiprocessing-pool-map-with-multiple-arguments
    repeated_di = [dataset_inputs] * num_fitted_replicas
    repeated_theoryid = [theoryid] * num_fitted_replicas
    repeated_fit = [str(fit)] * num_fitted_replicas
    repeated_mcseed = [mcseed] * num_fitted_replicas
    repeated_genrep = [genrep] * num_fitted_replicas
    # Import here to avoid cyclical imports, but we don't
    # want to import this in the task function because
    # loading the api is expensive.
    from validphys.api import API
    repeated_api = [API] * num_fitted_replicas

    args = zip(
        repeated_di,
        repeated_theoryid,
        repeated_fit,
        repeated_mcseed,
        repeated_genrep,
        repeated_api,
        replicas
    )

    if NPROC is None:
        NPROC = mp.cpu_count()
        log.warning(
            f"Using all {NPROC} cores available, this may be dangerous "
            "especially for use on a cluster. Consider setting the NPROC "
            "variable to something sensible."
        )

    if NPROC == 1:
        res = list(map(lambda x: _make_replica_task(*x), args))
    else:
        with mp.Pool(processes=NPROC) as pool:
            res = pool.starmap(_make_replica_task, args)

    df = pd.concat(res, axis=1)
    df.columns = [f"replica {i}" for i in replicas]

    return df


def recreate_pdf_pseudodata(pdf, NPROC=None):
    """Function that recreates the pseudodata of a PDF set
    accounting for the postfit reordering.

    # TODO: finish this example
    Example
    -------
    """
    # Import here to avoid cyclical imports
    from validphys.api import API

    fit_pseudodata = API.recreate_fit_pseudodata(fit=str(pdf), NPROC=NPROC)
    fitted_replica_indexes = API.fitted_replica_indexes(pdf=str(pdf))
    cols = [f"replica {i}" for i in fitted_replica_indexes]
    return fit_pseudodata.loc[:, cols]


def _datasets_mask(experiment_list):
    """Function to obtain a per datasets training/validation
    mask given the mask for the corresponding experiment.

    Returns
    -------
    dict:
        - tr_mask: training mask for the datasets in the experiment
        - vl_mask: validation mask for the datasets in the experiment
    """
    tr_mask = experiment_list["trmask"]
    vl_mask = experiment_list["vlmask"]
    slices = []
    start = 0
    for i in experiment_list["datasets"]:
        ndata = i["ndata"]
        slices.append(start + ndata)
        start += ndata

    return {
        "trmask": np.split(tr_mask, slices[:-1]),
        "vlmask": np.split(vl_mask, slices[:-1]),
    }


def training_validation_pseudodata(get_pseudodata):
    """Generator to yield a dictionary of training and validation DataFrame
    per replica indexed appropriately using a MultiIndex
    """
    exp_infos = get_pseudodata
    columns = ["experiment", "dataset", "id"]
    # Loop over all initial replicas
    for replica in exp_infos:
        tr_records, tr_central_values = [], []
        vl_records, vl_central_values = [], []
        # Loop over experiments in given replica
        for experiment in replica:
            split_masks = _datasets_mask(experiment)
            tr_mask, vl_mask = split_masks["trmask"], split_masks["vlmask"]
            # While we're here extend the central_values of the experiment
            tr_central_values.extend(np.squeeze(experiment["expdata"]))
            vl_central_values.extend(np.squeeze(experiment["expdata_vl"]))
            # Loop over datasets in experiment
            for i, dataset in enumerate(experiment["datasets"]):
                tr_dataset_mask = tr_mask[i]
                vl_dataset_mask = vl_mask[i]
                tr_indices = np.array((range(dataset["ndata"])))[tr_dataset_mask]
                vl_indices = np.array((range(dataset["ndata"])))[vl_dataset_mask]
                for tr_idat in tr_indices:
                    tr_records.append(
                        dict(
                            [
                                ("experiment", experiment["name"]),
                                ("dataset", dataset["name"]),
                                ("id", tr_idat),
                            ]
                        )
                    )
                for vl_idat in vl_indices:
                    vl_records.append(
                        dict(
                            [
                                ("experiment", experiment["name"]),
                                ("dataset", dataset["name"]),
                                ("id", vl_idat),
                            ]
                        )
                    )

        tr_df = pd.DataFrame(tr_records, columns=columns)
        vl_df = pd.DataFrame(vl_records, columns=columns)

        tr_df.set_index(columns, inplace=True)
        vl_df.set_index(columns, inplace=True)

        tr_index = tr_df.index
        vl_index = vl_df.index
        tr_vl_dict = {
            "trdata": pd.DataFrame(tr_central_values, index=tr_index, columns=["data"]),
            "vldata": pd.DataFrame(vl_central_values, index=vl_index, columns=["data"]),
        }
        yield tr_vl_dict
