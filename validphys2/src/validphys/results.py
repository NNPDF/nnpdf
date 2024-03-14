"""
results.py

Tools to obtain theory predictions and basic statistical estimators.
"""
from __future__ import generator_stop

from collections import OrderedDict, namedtuple
from collections.abc import Sequence
import logging

import numpy as np
import pandas as pd
import scipy.linalg as la

from reportengine import collect
from reportengine.checks import check_not_empty, remove_outer, require_one
from reportengine.table import table
from validphys.calcutils import all_chi2, bootstrap_values, calc_chi2, calc_phi, central_chi2
from validphys.checks import (
    check_pdf_is_montecarlo,
    check_speclabels_different,
    check_two_dataspecs,
)
from validphys.convolution import PredictionsRequireCutsError, central_predictions, predictions
from validphys.core import PDF, DataGroupSpec, DataSetSpec, Stats
from validphys.plotoptions.core import get_info

log = logging.getLogger(__name__)


class Result:
    pass


class StatsResult(Result):
    def __init__(self, stats):
        self.stats = stats

    @property
    def rawdata(self):
        """Returns the raw data with shape (Npoints, Npdf)"""
        return self.stats.data.T

    @property
    def error_members(self):
        """Returns the error members with shape (Npoints, Npdf)"""
        return self.stats.error_members().T

    @property
    def central_value(self):
        return self.stats.central_value()

    @property
    def std_error(self):
        return self.stats.std_error()

    def __len__(self):
        """Returns the number of data points in the result"""
        return self.rawdata.shape[0]


class DataResult(StatsResult):
    """Holds the relevant information from a given dataset"""

    def __init__(self, dataset, covmat, sqrtcovmat):
        loaded_cd = dataset.load_commondata()
        if isinstance(loaded_cd, list):
            cv = np.concatenate([cd.get_cv() for cd in loaded_cd])
        else:
            cv = loaded_cd.get_cv()
        self._central_value = cv
        stats = Stats(self._central_value)
        self._covmat = covmat
        self._sqrtcovmat = sqrtcovmat
        self._dataset = dataset
        super().__init__(stats)

    @property
    def label(self):
        return "Data"

    @property
    def central_value(self):
        return self._central_value

    @property
    def std_error(self):
        return np.sqrt(np.diag(self.covmat))

    @property
    def covmat(self):
        return self._covmat

    @property
    def sqrtcovmat(self):
        """Lower part of the Cholesky decomposition"""
        return self._sqrtcovmat

    @property
    def name(self):
        return self._dataset.name


class ThPredictionsResult(StatsResult):
    """Class holding theory prediction, inherits from StatsResult
    When created with `from_convolution`, it keeps tracks of the PDF for which it was computed
    """

    def __init__(
        self, dataobj, stats_class, datasetnames=None, label=None, pdf=None, theoryid=None
    ):
        self.stats_class = stats_class
        self.label = label
        self._datasetnames = datasetnames
        statsobj = stats_class(dataobj.T)
        self._pdf = pdf
        self._theoryid = theoryid
        super().__init__(statsobj)

    @staticmethod
    def make_label(pdf, dataset):
        """Deduce a reasonable label for the result based on pdf and dataspec"""
        th = dataset.thspec
        if hasattr(pdf, "label"):
            if hasattr(th, "label"):
                label = " ".join((pdf.label, th.label))
            else:
                label = pdf.label
        elif hasattr(th, "label"):
            label = th.label
        else:
            label = "{}@<Theory {}>".format(pdf, th.id)
        return label

    @classmethod
    def from_convolution(cls, pdf, dataset, central_only=False):
        # This should work for both single dataset and whole groups
        try:
            datasets = dataset.datasets
        except AttributeError:
            datasets = (dataset,)

        try:
            if central_only:
                preds = [central_predictions(d, pdf) for d in datasets]
            else:
                preds = [predictions(d, pdf) for d in datasets]
            th_predictions = pd.concat(preds)
        except PredictionsRequireCutsError as e:
            raise PredictionsRequireCutsError(
                "Predictions from FKTables always require cuts, "
                "if you want to use the fktable intrinsic cuts set `use_cuts: 'internal'`"
            ) from e

        label = cls.make_label(pdf, dataset)
        thid = dataset.thspec.id
        datasetnames = [i.name for i in datasets]
        return cls(th_predictions, pdf.stats_class, datasetnames, label, pdf=pdf, theoryid=thid)

    @property
    def datasetnames(self):
        return self._datasetnames


class ThUncertaintiesResult(StatsResult):
    """Class holding central theory predictions and the error bar corresponding to
    the theory uncertainties considered.
    The error members of this class correspond to central +- error_bar
    """

    def __init__(self, central, std_err, label=None):
        # All should be (ndata, 1)
        self._central = central
        self._std_err = std_err
        self.stats = None
        self.label = label

    @property
    def rawdata(self):
        return self._central

    @property
    def error_members(self):
        upper = self._central + self._std_err
        lower = self._central - self._std_err
        return np.concatenate([self._central, lower, upper], axis=1)

    @property
    def central_value(self):
        return self._central

    @property
    def std_error(self):
        return self._std_err

    def __len__(self):
        return self._central.shape[0]


class PositivityResult(StatsResult):
    @classmethod
    def from_convolution(cls, pdf, posset):
        data = predictions(posset, pdf)
        stats = pdf.stats_class(data.T)
        return cls(stats)


def data_index(data):
    """
    Given a core.DataGroupSpec instance, return pd.MultiIndex
    with the following levels:

    1. experiment
    2. datasets
    3. datapoints indices (cuts already applied to)


    Parameters
    ----------
    data: core.DataGroupSpec

    Returns
    -------
    pd.MultiIndex

    """
    tuples = []
    for ds in data.datasets:
        experiment = get_info(ds).experiment

        for i in ds.cuts.load():
            tp = (experiment, ds.name, i)
            tuples.append(tp)
    return pd.MultiIndex.from_tuples(tuples, names=('experiment', 'dataset', 'id'))


# TODO: finish deprecating all dependencies on this index largely in theorycovmat module
groups_data = collect("data", ("group_dataset_inputs_by_metadata",))

experiments_data = collect("data", ("group_dataset_inputs_by_experiment",))

procs_data = collect("data", ("group_dataset_inputs_by_process",))


def groups_index(groups_data):
    """Return a pandas.MultiIndex with levels for group, dataset and point
    respectively, the group is determined by a key in the dataset metadata, and
    controlled by `metadata_group` key in the runcard.

    Example
    -------
    TODO: add example

    """
    records = []
    for group in groups_data:
        for dataset in group.datasets:
            if dataset.cuts:
                data_id = dataset.cuts.load()
            else:
                # No cuts - use all data
                data_id = np.arange(dataset.commondata.ndata, dtype=int)
            for idat in data_id:
                records.append(
                    dict([("group", str(group.name)), ("dataset", str(dataset.name)), ("id", idat)])
                )

    columns = ["group", "dataset", "id"]
    df = pd.DataFrame(records, columns=columns)
    df.set_index(columns, inplace=True)
    return df.index


def experiments_index(experiments_data):
    return groups_index(experiments_data)


def procs_index(procs_data):
    return groups_index(procs_data)


def groups_data_values(group_result_table):
    """Returns list of data values for the input groups."""
    data_central_values = group_result_table["data_central"]
    return data_central_values


def procs_data_values(proc_result_table):
    """Like groups_data_values but grouped by process."""
    data_central_values = proc_result_table["data_central"]
    return data_central_values


def procs_data_values_experiment(proc_result_table_experiment):
    """Like groups_data_values but grouped by experiment."""
    data_central_values = proc_result_table_experiment["data_central"]
    return data_central_values


groups_results = collect("dataset_inputs_results", ("group_dataset_inputs_by_metadata",))

procs_results = collect("dataset_inputs_results_central", ("group_dataset_inputs_by_process",))

procs_results_experiment = collect(
    "dataset_inputs_results_central", ("group_dataset_inputs_by_experiment",)
)

groups_results_central = collect(
    "dataset_inputs_results_central", ("group_dataset_inputs_by_metadata",)
)


def group_result_central_table_no_table(groups_results_central, groups_index):
    """Generate a table containing the data central value and the central prediction"""
    result_records = []
    for group_results in groups_results_central:
        dt, th = group_results
        for dt_central, th_central in zip(dt.central_value, th.central_value):
            result_records.append(
                dict([("data_central", dt_central), ("theory_central", th_central)])
            )
    if not result_records:
        log.warning("Empty records for group results")
        return pd.DataFrame()
    df = pd.DataFrame(result_records, columns=result_records[0].keys(), index=groups_index)

    return df


def group_result_table_no_table(groups_results, groups_index):
    """Generate a table containing the data central value, the central prediction,
    and the prediction for each PDF member."""
    result_records = []
    for group_results in groups_results:
        dt, th = group_results
        for index, (dt_central, th_central) in enumerate(zip(dt.central_value, th.central_value)):
            replicas = (
                ("rep_%05d" % (i + 1), th_rep)
                for i, th_rep in enumerate(th.error_members[index, :])
            )

            result_records.append(
                dict([("data_central", dt_central), ("theory_central", th_central), *replicas])
            )
    if not result_records:
        log.warning("Empty records for group results")
        return pd.DataFrame()
    df = pd.DataFrame(result_records, columns=result_records[0].keys(), index=groups_index)

    return df


@table
def group_result_table(group_result_table_no_table):
    """Duplicate of group_result_table_no_table but with a table decorator."""
    return group_result_table_no_table


def proc_result_table_no_table(procs_results, procs_index):
    return group_result_table_no_table(procs_results, procs_index)


@table
def proc_result_table(proc_result_table_no_table):
    return proc_result_table_no_table


@table
def proc_result_table_experiment(procs_results_experiment, experiments_index):
    return group_result_table_no_table(procs_results_experiment, experiments_index)


experiment_result_table = collect("group_result_table", ("group_dataset_inputs_by_experiment",))


@table
def group_result_table_68cl(groups_results, group_result_table_no_table: pd.DataFrame, pdf: PDF):
    """Generate a table containing the data central value, the data 68% confidence levels, the central prediction,
    and 68% confidence level bounds of the prediction.
    """
    df = group_result_table_no_table
    # replica data is every columns after central values, transpose for stats class
    replica_data = df.iloc[:, 2:].values.T
    # Use pdf stats class but reshape output to have each row as a data point
    th_unc_array = [level.reshape(-1, 1) for level in pdf.stats_class(replica_data).errorbar68()]
    # concatenate for dataframe construction
    th_unc_array_reshaped = np.concatenate(th_unc_array, axis=1)
    data_unc_array = np.concatenate([i[0].std_error for i in groups_results])
    uncertainties_array = np.c_[data_unc_array, th_unc_array_reshaped]
    df_cl = pd.DataFrame(
        uncertainties_array,
        index=df.index,
        columns=["data uncertainty", "theory_lower", "theory_upper"],
    )
    res = pd.concat([df.iloc[:, :2], df_cl], axis=1)
    return res


experiments_covmat_collection = collect(
    "dataset_inputs_covariance_matrix", ("group_dataset_inputs_by_experiment",)
)


def experiments_covmat_no_table(experiments_data, experiments_index, experiments_covmat_collection):
    """Makes the total experiments covariance matrix, which can then
    be reindexed appropriately by the chosen grouping. The covariance
    matrix must first be grouped by experiments to ensure correlations
    within experiments are preserved."""
    data = np.zeros((len(experiments_index), len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, experiment_covmat in zip(experiments_data, experiments_covmat_collection):
        name = experiment.name
        df.loc[[name], [name]] = experiment_covmat
    return df


def relabel_experiments_to_groups(input_covmat, groups_index):
    """Takes a covmat grouped by experiments and relabels
    it by groups. This allows grouping over experiments to
    preserve experimental correlations outwith the chosen
    grouping."""
    # Sorting along dataset axis so we can apply the groups index directly
    input_covmat = input_covmat.sort_index(axis=0, level=1)
    input_covmat = input_covmat.sort_index(axis=1, level=1)
    sorted_groups_index = groups_index.sortlevel(1)[0]
    df = pd.DataFrame(input_covmat.values, index=sorted_groups_index, columns=sorted_groups_index)
    # Reindexing to fit with groups_index
    df = df.reindex(groups_index, axis=0)
    df = df.reindex(groups_index, axis=1)
    return df


def groups_covmat_no_table(experiments_covmat_no_table, groups_index):
    """Export the covariance matrix for the groups. It exports the full
    (symmetric) matrix, with the 3 first rows and columns being:

        - group name

        - dataset name

        - index of the point within the dataset.
    """
    return relabel_experiments_to_groups(experiments_covmat_no_table, groups_index)


@table
def groups_covmat(groups_covmat_no_table):
    """Duplicate of groups_covmat_no_table but with a table decorator."""
    return groups_covmat_no_table


def procs_covmat_no_table(experiments_covmat_no_table, procs_index):
    return relabel_experiments_to_groups(experiments_covmat_no_table, procs_index)


@table
def procs_covmat(procs_covmat_no_table):
    return procs_covmat_no_table


experiments_sqrt_covmat = collect(
    "dataset_inputs_sqrt_covmat", ("group_dataset_inputs_by_experiment",)
)


@table
def experiments_sqrtcovmat(experiments_data, experiments_index, experiments_sqrt_covmat):
    """Like experiments_covmat, but dump the lower triangular part of the
    Cholesky decomposition as used in the fit. The upper part indices are set
    to zero.
    """
    data = np.zeros((len(experiments_index), len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, experiments_sqrt_covmat in zip(experiments_data, experiments_sqrt_covmat):
        name = experiment.name
        experiments_sqrt_covmat[np.triu_indices_from(experiments_sqrt_covmat, k=1)] = 0
        df.loc[[name], [name]] = experiments_sqrt_covmat
    return df


@table
def groups_sqrtcovmat(experiments_sqrtcovmat, groups_index):
    """Like experiments_sqrtcovmat but relabelled to the chosen grouping."""
    return relabel_experiments_to_groups(experiments_sqrtcovmat, groups_index)


@table
def experiments_invcovmat(experiments_data, experiments_index, experiments_covmat_collection):
    """Compute and export the inverse covariance matrix.
    Note that this inverts the matrices with the LU method which is
    suboptimal."""
    data = np.zeros((len(experiments_index), len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, experiment_covmat in zip(experiments_data, experiments_covmat_collection):
        name = experiment.name
        # Improve this inversion if this method tuns out to be important
        invcov = la.inv(experiment_covmat)
        df.loc[[name], [name]] = invcov
    return df


@table
def groups_invcovmat(experiments_invcovmat, groups_index):
    """Like experiments_invcovmat but relabelled to the chosen grouping."""
    return relabel_experiments_to_groups(experiments_invcovmat, groups_index)


@table
def groups_normcovmat(groups_covmat, groups_data_values):
    """Calculates the grouped experimental covariance matrix normalised to data."""
    df = groups_covmat
    index = df.index
    # Reindexing data so that it is aligned with the covmat
    groups_data_values = groups_data_values.reindex(index)
    groups_data_array = np.array(groups_data_values)
    mat = df / np.outer(groups_data_array, groups_data_array)
    return mat


@table
def procs_normcovmat(procs_covmat, procs_data_values):
    return groups_normcovmat(procs_covmat, procs_data_values)


@table
def groups_corrmat(groups_covmat):
    """Generates the grouped experimental correlation matrix with groups_covmat as input"""
    df = groups_covmat
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat)) ** (-0.5)
    mat = diag_minus_half[:, np.newaxis] * df * diag_minus_half
    return mat


@table
def procs_corrmat(procs_covmat):
    return groups_corrmat(procs_covmat)


def results(dataset: (DataSetSpec), pdf: PDF, covariance_matrix, sqrt_covmat):
    """Tuple of data and theory results for a single pdf. The data will have an associated
    covariance matrix, which can include a contribution from the theory covariance matrix which
    is constructed from scale variation. The inclusion of this covariance matrix by default is used
    where available, however this behaviour can be modified with the flag `use_theorycovmat`.

    The theory is specified as part of the dataset (a remnant of the old C++ layout)
    A group of datasets is also allowed.
    """
    # TODO: is the message about the usage of the theory covariance matrix here true?
    # probably not in most cases...
    return (
        DataResult(dataset, covariance_matrix, sqrt_covmat),
        ThPredictionsResult.from_convolution(pdf, dataset),
    )


def results_central(dataset: (DataSetSpec), pdf: PDF, covariance_matrix, sqrt_covmat):
    """Same as :py:func:`results` but only calculates the prediction for replica0."""
    return (
        DataResult(dataset, covariance_matrix, sqrt_covmat),
        ThPredictionsResult.from_convolution(pdf, dataset, central_only=True),
    )


def results_with_theory_covmat(dataset, results, theory_covmat_dataset):
    """Returns results with a modfy ``DataResult`` such that the covariance matrix includes
    also the theory covmat.
    This can be used to make use of results that consider scale variations without including
    the theory covmat as part of the covariance matrix used by other validphys function.
    Most notably, this can be used to compute the chi2 including theory errors while plotting
    data theory covariance in which the experimental uncertainties are not stained by the thcovmat
    """
    # TODO: in principle this function could be removed, and `results` could automagically
    # include the theory covmat when `use_theorycovmat: true` by changing the nodes in `config.py`
    # however at the moment config.py _loads_ theory covmats and we need to compute it on the fly
    from .covmats import sqrt_covmat

    data_result, central_th_result = results
    total_covmat = theory_covmat_dataset + data_result.covmat

    data_result = DataResult(dataset, total_covmat, sqrt_covmat(total_covmat))
    return (data_result, central_th_result)


def results_with_scale_variations(results, theory_covmat_dataset):
    """Use the theory covariance matrix to generate a ThPredictionsResult-compatible object
    modified so that its uncertainties correspond to a combination of
    the PDF and theory (scale variations) errors added in quadrature.
    This allows to plot results including scale variations

    By doing this we lose all information about prediction for the individual replicas or theories
    """
    data_result, central_th_result = results

    # Use the central value and PDF error from the central theory
    cv = central_th_result.central_value
    pdf_error = central_th_result.std_error
    sv_error_sq = np.diag(theory_covmat_dataset)

    total_error = np.sqrt(pdf_error**2 + sv_error_sq)

    theory_error_result = ThUncertaintiesResult(cv, total_error, label=central_th_result.label)
    return (data_result, theory_error_result)


def dataset_inputs_results_central(
    data, pdf: PDF, dataset_inputs_covariance_matrix, dataset_inputs_sqrt_covmat
):
    """Like `dataset_inputs_results` but for a group of datasets and replica0."""
    return results_central(data, pdf, dataset_inputs_covariance_matrix, dataset_inputs_sqrt_covmat)


def dataset_inputs_results(
    data, pdf: PDF, dataset_inputs_covariance_matrix, dataset_inputs_sqrt_covmat
):
    """Like `results` but for a group of datasets"""
    return results(data, pdf, dataset_inputs_covariance_matrix, dataset_inputs_sqrt_covmat)


# It's better to duplicate a few lines than to complicate the logic of
# ``results`` to support this.
# TODO: The above comment doesn't make sense after adding T0. Deprecate this
def pdf_results(
    dataset: (DataSetSpec, DataGroupSpec), pdfs: Sequence, covariance_matrix, sqrt_covmat
):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    th_results = [ThPredictionsResult.from_convolution(pdf, dataset) for pdf in pdfs]

    return (DataResult(dataset, covariance_matrix, sqrt_covmat), *th_results)


@require_one("pdfs", "pdf")
@remove_outer("pdfs", "pdf")
def one_or_more_results(
    dataset: (DataSetSpec, DataGroupSpec),
    covariance_matrix,
    sqrt_covmat,
    pdfs: (type(None), Sequence) = None,
    pdf: (type(None), PDF) = None,
):
    """Generate a list of results, where the first element is the data values,
    and the next is either the prediction for pdf or for each of the pdfs.
    Which of the two is selected intelligently depending on the namespace,
    when executing as an action."""
    if pdf is not None:
        return results(dataset, pdf, covariance_matrix, sqrt_covmat)
    return pdf_results(dataset, pdfs, covariance_matrix, sqrt_covmat)


Chi2Data = namedtuple("Chi2Data", ("replica_result", "central_result", "ndata"))


def abs_chi2_data(results):
    """Return a tuple (member_chi², central_chi², numpoints) for a
    given dataset"""
    data_result, th_result = results

    chi2s = all_chi2(results)

    central_result = central_chi2(results)

    return Chi2Data(th_result.stats_class(chi2s[:, np.newaxis]), central_result, len(data_result))


def abs_chi2_data_thcovmat(results_with_theory_covmat):
    """The same as ``abs_chi2_data`` but considering as well the theory uncertainties"""
    return abs_chi2_data(results_with_theory_covmat)


def dataset_inputs_abs_chi2_data(dataset_inputs_results):
    """Like `abs_chi2_data` but for a group of inputs"""
    return abs_chi2_data(dataset_inputs_results)


def phi_data(abs_chi2_data):
    """Calculate phi using values returned by `abs_chi2_data`.

    Returns tuple of (float, int): (phi, numpoints)

    For more information on how phi is calculated see Eq.(24) in
    1410.8849
    """
    alldata, central, npoints = abs_chi2_data
    return (np.sqrt((alldata.error_members().mean() - central) / npoints), npoints)


def dataset_inputs_phi_data(dataset_inputs_abs_chi2_data):
    """Like `phi_data` but for group of datasets"""
    return phi_data(dataset_inputs_abs_chi2_data)


experiments_phi_data = collect("dataset_inputs_phi_data", ("group_dataset_inputs_by_experiment",))


def total_phi_data_from_experiments(experiments_phi_data):
    """Like :py:func:`dataset_inputs_phi_data` except calculate phi for
    each experiment and then sum the contributions. Note that since
    the definition of phi is

        phi = sqrt( (<chi2[T_k]> - chi2[<T_k>]) / n_data ),

    where k is the replica index, the total phi is

        sqrt( sum(n_data*phi**2) / sum(n_data) )

    where the sums run over experiment

    This is only a valid method of calculating total phi provided that there are
    no inter-experimental correlations.

    """

    unnorm_phi_squared, ndata = np.sum(
        [(ndata * phi**2, ndata) for phi, ndata in experiments_phi_data], axis=0
    )
    return np.sqrt(unnorm_phi_squared / ndata), ndata


@check_pdf_is_montecarlo
def dataset_inputs_bootstrap_phi_data(dataset_inputs_results, bootstrap_samples=500):
    """Takes the data result and theory prediction given `dataset_inputs` and
    then returns a bootstrap distribution of phi.
    By default `bootstrap_samples` is set to a sensible value (500). However
    a different value can be specified in the runcard.

    For more information on how phi is calculated see `phi_data`
    """
    dt, th = dataset_inputs_results
    diff = np.array(th.error_members - dt.central_value[:, np.newaxis])
    phi_resample = bootstrap_values(
        diff, bootstrap_samples, apply_func=(lambda x, y: calc_phi(y, x)), args=[dt.sqrtcovmat]
    )
    return phi_resample


@check_pdf_is_montecarlo
def dataset_inputs_bootstrap_chi2_central(
    dataset_inputs_results, bootstrap_samples=500, boot_seed=123
):
    """Takes the data result and theory prediction given dataset_inputs and
    then returns a bootstrap distribution of central chi2.
    By default `bootstrap_samples` is set to a sensible value (500). However
    a different value can be specified in the runcard.
    """
    dt, th = dataset_inputs_results
    diff = np.array(th.error_members - dt.central_value[:, np.newaxis])
    cchi2 = lambda x, y: calc_chi2(y, x.mean(axis=1))
    chi2_central_resample = bootstrap_values(
        diff, bootstrap_samples, boot_seed=boot_seed, apply_func=(cchi2), args=[dt.sqrtcovmat]
    )
    return chi2_central_resample


@table
def predictions_by_kinematics_table(results, kinematics_table_notable):
    """Return a table combining the output of
    :py:func:`validphys.kinematics.kinematics_table`` with the data and theory
    central values."""
    tb = kinematics_table_notable.copy()
    data, theory = results
    tb['data'] = data.central_value
    tb['prediction'] = theory.central_value
    return tb


groups_each_dataset_chi2 = collect("each_dataset_chi2", ("group_dataset_inputs_by_metadata",))
groups_chi2_by_process = collect(
    "dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_process",)
)
groups_each_dataset_chi2_by_process = collect(
    "each_dataset_chi2", ("group_dataset_inputs_by_process",)
)


@table
def groups_chi2_table(groups_data, pdf, groups_chi2, groups_each_dataset_chi2):
    """Return a table with the chi² to the groups and each dataset in
    the groups, grouped by metadata."""
    records = []
    for group, groupres, dsresults in zip(groups_data, groups_chi2, groups_each_dataset_chi2):
        for dataset, dsres in zip(group, dsresults):
            stats = chi2_stats(dsres)
            stats["group"] = dataset.name
            records.append(stats)
    return pd.DataFrame(records)


experiments_chi2_table = collect("groups_chi2_table", ("group_dataset_inputs_by_experiment",))


@table
def procs_chi2_table(procs_data, pdf, groups_chi2_by_process, groups_each_dataset_chi2_by_process):
    """Same as groups_chi2_table but by process"""
    return groups_chi2_table(
        procs_data, pdf, groups_chi2_by_process, groups_each_dataset_chi2_by_process
    )


def positivity_predictions_data_result(pdf, posdataset):
    """Return an object containing the values of the positivuty observable."""
    return PositivityResult.from_convolution(pdf, posdataset)


positivity_predictions_for_pdfs = collect(positivity_predictions_data_result, ("pdfs",))
dataspecs_positivity_predictions = collect(positivity_predictions_data_result, ("dataspecs",))
dataspecs_posdataset = collect("posdataset", ("dataspecs",))


def count_negative_points(possets_predictions):
    """Return the number of replicas with negative predictions for each bin
    in the positivity observable."""
    return np.sum([(r.error_members < 0).sum(axis=0) for r in possets_predictions], axis=0)


chi2_stat_labels = {
    "central_mean": r"$\chi^2_{rep0}$",
    "npoints": r"$N_{data}$",
    "perreplica_mean": r"$\left< \chi^2 \right>_{rep}$",
    "perreplica_std": r"$std_{rep}(\chi^2)$",
    "chi2_per_data": r"$\chi^2 / N_{data}$",
}


def experiments_chi2_stats(total_chi2_data):
    """Compute several estimators from the chi² for an
    aggregate of experiments:

     - central_mean

     - npoints

     - perreplica_mean

     - perreplica_std

     - chi2_per_data
    """
    rep_data, central_result, npoints = total_chi2_data
    m = central_result.mean()
    rep_mean = rep_data.error_members().mean()
    return OrderedDict(
        [
            ("central_mean", m),
            ("npoints", npoints),
            ("chi2_per_data", m / npoints),
            ("perreplica_mean", rep_mean),
            ("perreplica_std", rep_data.std_error().mean()),
        ]
    )


def chi2_stats(abs_chi2_data):
    """Compute several estimators from the chi²:

    - central_mean

    - npoints

    - perreplica_mean

    - perreplica_std

    - chi2_per_data
    """
    rep_data, central_result, npoints = abs_chi2_data
    m = central_result.mean()
    rep_mean = rep_data.error_members().mean()
    return OrderedDict(
        [
            ("central_mean", m),
            ("npoints", npoints),
            ("chi2_per_data", m / npoints),
            ("perreplica_mean", rep_mean),
            ("perreplica_std", rep_data.std_error().mean()),
        ]
    )


@table
def dataset_chi2_table(chi2_stats, dataset):
    """Show the chi² estimators for a given dataset"""
    return pd.DataFrame(chi2_stats, index=[dataset.name])


groups_chi2 = collect("dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_metadata",))

procs_chi2 = collect("dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_process",))

fits_groups_chi2_data = collect("groups_chi2", ("fits", "fitcontext"))
fits_groups = collect("groups_data", ("fits", "fitcontext"))


# TODO: Possibly get rid of the per_point_data parameter and have separate
# actions for absolute and relative tables.
@table
def fits_groups_chi2_table(
    fits_name_with_covmat_label, fits_groups, fits_groups_chi2_data, per_point_data: bool = True
):
    """A table with the chi2 computed with the theory corresponding to each fit
    for all datasets in the fit, grouped according to a key in the metadata, the
    grouping can be controlled with `metadata_group`.

    If points_per_data is True, the chi² will be shown divided by ndata.
    Otherwise chi² values will be absolute.

    """
    dfs = []
    cols = ("ndata", r"$\chi^2/ndata$") if per_point_data else ("ndata", r"$\chi^2$")
    for label, groups, groups_chi2 in zip(
        fits_name_with_covmat_label, fits_groups, fits_groups_chi2_data
    ):
        records = []
        for group, group_chi2 in zip(groups, groups_chi2):
            mean_chi2 = group_chi2.central_result.mean()
            npoints = group_chi2.ndata
            records.append(dict(group=str(group), npoints=npoints, mean_chi2=mean_chi2))
        df = pd.DataFrame.from_records(
            records, columns=("group", "npoints", "mean_chi2"), index=("group",)
        )
        if per_point_data:
            df["mean_chi2"] /= df["npoints"]
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    res = pd.concat(dfs, axis=1)
    return res


groups_phi = collect("dataset_inputs_phi_data", ("group_dataset_inputs_by_metadata",))
fits_groups_phi = collect("groups_phi", ("fits", "fitcontext"))


@table
def fits_groups_phi_table(fits_name_with_covmat_label, fits_groups, fits_groups_phi):
    """For every fit, returns phi and number of data points for each group of
    datasets, which are grouped according to a key in the metadata. The behaviour
    of the grouping can be controlled with `metadata_group` runcard key.

    """
    dfs = []
    cols = ("ndata", r"$\phi$")
    for label, groups, groups_phi in zip(fits_name_with_covmat_label, fits_groups, fits_groups_phi):
        records = []
        for group, (group_phi, npoints) in zip(groups, groups_phi):
            records.append(dict(group=str(group), npoints=npoints, phi=group_phi))
        df = pd.DataFrame.from_records(
            records, columns=("group", "npoints", "phi"), index=("group",)
        )
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    res = pd.concat(dfs, axis=1)
    return res


@table
@check_speclabels_different
def dataspecs_groups_chi2_table(
    dataspecs_speclabel, dataspecs_groups, dataspecs_groups_chi2_data, per_point_data: bool = True
):
    """Same as fits_groups_chi2_table but for an arbitrary list of dataspecs."""
    return fits_groups_chi2_table(
        dataspecs_speclabel,
        dataspecs_groups,
        dataspecs_groups_chi2_data,
        per_point_data=per_point_data,
    )


# we need this to reorder the datasets correctly, a potentially more satisfactory
# method could be to make a datasets chi2 table which gets collected and concatenated
groups_datasets_chi2_data = collect("each_dataset_chi2", ("group_dataset_inputs_by_metadata",))
fits_datasets_chi2_data = collect("groups_datasets_chi2_data", ("fits", "fitcontext"))


@table
def fits_datasets_chi2_table(
    fits_name_with_covmat_label, fits_groups, fits_datasets_chi2_data, per_point_data: bool = True
):
    """A table with the chi2 for each included dataset in the fits, computed
    with the theory corresponding to the fit. The result are indexed in two
    levels by experiment and dataset, where experiment is the grouping of datasets according to the
    `experiment` key in the PLOTTING info file.  If points_per_data is True, the chi² will be shown
    divided by ndata. Otherwise they will be absolute."""

    cols = ("ndata", r"$\chi^2/ndata$") if per_point_data else ("ndata", r"$\chi^2$")

    dfs = []
    for label, groups, groups_dsets_chi2 in zip(
        fits_name_with_covmat_label, fits_groups, fits_datasets_chi2_data
    ):
        records = []
        for group, dsets_chi2 in zip(groups, groups_dsets_chi2):
            for dataset, chi2 in zip(group.datasets, dsets_chi2):
                ndata = chi2.ndata

                records.append(
                    dict(
                        group=str(group),
                        dataset=str(dataset),
                        npoints=ndata,
                        mean_chi2=chi2.central_result.mean(),
                    )
                )

        df = pd.DataFrame.from_records(
            records,
            columns=("group", "dataset", "npoints", "mean_chi2"),
            index=("group", "dataset"),
        )
        if per_point_data:
            df["mean_chi2"] /= df["npoints"]
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    return pd.concat(dfs, axis=1)


dataspecs_datasets_chi2_data = collect("groups_datasets_chi2_data", ("dataspecs",))


@table
@check_speclabels_different
def dataspecs_datasets_chi2_table(
    dataspecs_speclabel, dataspecs_groups, dataspecs_datasets_chi2_data, per_point_data: bool = True
):
    """Same as fits_datasets_chi2_table but for arbitrary dataspecs."""
    return fits_datasets_chi2_table(
        dataspecs_speclabel,
        dataspecs_groups,
        dataspecs_datasets_chi2_data,
        per_point_data=per_point_data,
    )


fits_total_chi2_data = collect("total_chi2_data", ("fits", "fitcontext"))
dataspecs_total_chi2_data = collect("total_chi2_data", ("dataspecs",))


# TODO: Decide what to do with the horrible totals code.
@table
def fits_chi2_table(
    fits_total_chi2_data, fits_datasets_chi2_table, fits_groups_chi2_table, show_total: bool = False
):
    """Show the chi² of each and number of points of each dataset and experiment
    of each fit, where experiment is a group of datasets according to the `experiment` key in
    the PLOTTING info file, computed with the theory corresponding to the fit. Dataset that are not
    included in some fit appear as `NaN`
    """
    lvs = fits_groups_chi2_table.index
    # The explicit call to list is because pandas gets confused otherwise
    expanded_index = pd.MultiIndex.from_product((list(lvs), ["Total"]))
    edf = fits_groups_chi2_table.set_index(expanded_index)
    ddf = fits_datasets_chi2_table
    dfs = []
    # TODO: Better way to do the merge preserving the order?
    for lv in lvs:
        dfs.append(pd.concat((edf.loc[lv], ddf.loc[lv]), copy=False, axis=0))
    if show_total:
        total_points = np.array([total_chi2_data.ndata for total_chi2_data in fits_total_chi2_data])
        total_chi = np.array(
            [total_chi2_data.central_result for total_chi2_data in fits_total_chi2_data]
        )
        total_chi /= total_points
        row = np.zeros(len(total_points) * 2)
        row[::2] = total_points
        row[1::2] = total_chi
        df = pd.DataFrame(
            np.atleast_2d(row), columns=fits_groups_chi2_table.columns, index=["Total"]
        )
        dfs.append(df)
        keys = [*lvs, "Total"]
    else:
        keys = lvs

    res = pd.concat(dfs, axis=0, keys=keys)
    return res


@table
def dataspecs_chi2_table(
    dataspecs_total_chi2_data,
    dataspecs_datasets_chi2_table,
    dataspecs_groups_chi2_table,
    show_total: bool = False,
):
    """Same as fits_chi2_table but for an arbitrary list of dataspecs"""
    return fits_chi2_table(
        dataspecs_total_chi2_data,
        dataspecs_datasets_chi2_table,
        dataspecs_groups_chi2_table,
        show_total,
    )


@table
@check_two_dataspecs
def dataspecs_chi2_differences_table(dataspecs, dataspecs_chi2_table):
    """Given two dataspecs, print the chi² (using dataspecs_chi2_table)
    and the difference between the first and the second."""
    df = dataspecs_chi2_table.copy()
    # TODO: Make this mind the number of points somehow
    diff = df.iloc[:, 1] - df.iloc[:, 3]
    df["difference"] = diff
    return df


experiments_chi2_data = collect(
    "dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_experiment",)
)


def total_chi2_data_from_experiments(experiments_chi2_data, pdf):
    """Like :py:func:`dataset_inputs_abs_chi2_data`, except sums the contribution
    from each experiment which is more efficient in the case that the total
    covariance matrix is block diagonal in experiments.

    This is valid as long as there are no cross experiment correlations from
    e.g. theory covariance matrices.

    """

    central_result = np.sum(
        [exp_chi2_data.central_result for exp_chi2_data in experiments_chi2_data]
    )

    # we sum data, not error_members here because we feed it back into the stats
    # class, the stats class error_members cuts off the CV if needed
    data_sum = np.sum(
        [exp_chi2_data.replica_result.data for exp_chi2_data in experiments_chi2_data], axis=0
    )

    ndata = np.sum([exp_chi2_data.ndata for exp_chi2_data in experiments_chi2_data])
    return Chi2Data(pdf.stats_class(data_sum), central_result, ndata)


def dataset_inputs_chi2_per_point_data(dataset_inputs_abs_chi2_data):
    """Return the total chi²/ndata for all data, specified by dataset_inputs.
    Covariance matrix is fully correlated across datasets, with all known
    correlations.
    """
    return dataset_inputs_abs_chi2_data.central_result / dataset_inputs_abs_chi2_data.ndata


def total_chi2_per_point_data(total_chi2_data):
    return dataset_inputs_chi2_per_point_data(total_chi2_data)


@table
@check_not_empty("groups_data")
def perreplica_chi2_table(groups_data, groups_chi2, total_chi2_data):
    """Chi² per point for each replica for each group.
    Also outputs the total chi² per replica.
    The columns come in two levels: The first is the name of the group,
    and the second is the number of points."""

    chs = groups_chi2
    total_chis = np.zeros((len(groups_data) + 1, 1 + len(chs[0].replica_result.error_members())))
    ls = []
    for i, ch in enumerate(chs, 1):
        th, central, l = ch
        total_chis[i] = np.concatenate([[central], *th.error_members()])
        ls.append(l)
    total_rep, total_central, total_n = total_chi2_data
    total_chis[0] = np.concatenate([[total_central], *total_rep.error_members()])
    total_chis[0] /= total_n
    total_chis[1:, :] /= np.array(ls)[:, np.newaxis]

    columns = pd.MultiIndex.from_arrays(
        (["Total", *[str(exp) for exp in groups_data]], [total_n, *ls]), names=["name", "npoints"]
    )
    return pd.DataFrame(total_chis.T, columns=columns)


@table
def theory_description(theoryid):
    """A table with the theory settings."""
    return pd.DataFrame(pd.Series(theoryid.get_description()), columns=[theoryid])


def groups_central_values_no_table(group_result_central_table_no_table):
    """Returns a theoryid-dependent list of central theory predictions
    for a given group."""
    central_theory_values = group_result_central_table_no_table["theory_central"]
    return central_theory_values


@table
def groups_central_values(group_result_central_table_no_table):
    """Duplicate of groups_central_values_no_table but takes
    group_result_table rather than groups_central_values_no_table,
    and has a table decorator."""
    central_theory_values = group_result_central_table_no_table["theory_central"]
    return central_theory_values


def procs_central_values_no_table(proc_result_table_no_table):
    central_theory_values = proc_result_table_no_table["theory_central"]
    return central_theory_values


@table
def procs_central_values(procs_central_values_no_table):
    return procs_central_values_no_table


dataspecs_each_dataset_chi2 = collect("each_dataset_chi2", ("dataspecs",))
each_dataset = collect("dataset", ("data",))
dataspecs_each_dataset = collect("each_dataset", ("dataspecs",))


@table
@check_speclabels_different
def dataspecs_dataset_chi2_difference_table(
    dataspecs_each_dataset, dataspecs_each_dataset_chi2, dataspecs_speclabel
):
    r"""Returns a table with difference between the chi2 and the expected chi2
    in units of the expected chi2 standard deviation, given by

        chi2_diff = (\chi2 - N)/sqrt(2N)

    for each dataset for each dataspec.

    """
    dfs = []
    cols = [r"$(\chi^2 - N)/\sqrt{2N}$"]
    for label, datasets, chi2s in zip(
        dataspecs_speclabel, dataspecs_each_dataset, dataspecs_each_dataset_chi2
    ):
        records = []
        for dataset, chi2 in zip(datasets, chi2s):
            ndata = chi2.ndata

            records.append(
                dict(
                    dataset=str(dataset),
                    chi2_stat=(chi2.central_result.mean() - ndata) / np.sqrt(2 * ndata),
                )
            )

        df = pd.DataFrame.from_records(
            records, columns=("dataset", "chi2_stat"), index=("dataset",)
        )
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    return pd.concat(dfs, axis=1)


each_dataset_chi2 = collect(abs_chi2_data, ("data",))

pdfs_total_chi2 = collect(total_chi2_per_point_data, ("pdfs",))

# These are convenient ways to iterate and extract various data from fits
fits_chi2_data = collect(abs_chi2_data, ("fits", "fitcontext", "dataset_inputs"))

fits_total_chi2 = collect("total_chi2_per_point_data", ("fits", "fitcontext"))

fits_total_chi2_for_groups = collect("total_chi2_per_point_data", ("fits", "fittheoryandpdf"))

fits_pdf = collect("pdf", ("fits", "fitpdf"))


groups_data_phi = collect(dataset_inputs_phi_data, ("group_dataset_inputs_by_metadata",))
fits_groups_data_phi = collect("groups_data_phi", ("fits", "fitcontext"))
groups_bootstrap_phi = collect(
    dataset_inputs_bootstrap_phi_data, ("group_dataset_inputs_by_metadata",)
)

fits_groups_chi2 = collect("groups_chi2", ("fits", "fitcontext"))
fits_groups_data = collect("groups_data", ("fits", "fitcontext"))

# Collections over dataspecs
dataspecs_groups = collect("groups_data", ("dataspecs",))
dataspecs_groups_chi2_data = collect("groups_chi2", ("dataspecs",))
dataspecs_groups_bootstrap_phi = collect("groups_bootstrap_phi", ("dataspecs",))
dataspecs_results = collect("results", ("dataspecs",))
dataspecs_results_with_scale_variations = collect("results_with_scale_variations", ("dataspecs",))
dataspecs_total_chi2 = collect("total_chi2_per_point_data", ("dataspecs",))

dataspecs_speclabel = collect("speclabel", ("dataspecs",), element_default=None)
dataspecs_cuts = collect("cuts", ("dataspecs",))
dataspecs_experiments = collect("experiments", ("dataspecs",))
dataspecs_dataset = collect("dataset", ("dataspecs",))
dataspecs_commondata = collect("commondata", ("dataspecs",))
dataspecs_pdf = collect("pdf", ("dataspecs",))
dataspecs_fit = collect("fit", ("dataspecs",))
