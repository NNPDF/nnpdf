# -*- coding: utf-8 -*-
"""
results.py

Tools to obtain theory predictions and basic statistical estimators.
"""
from __future__ import generator_stop

from collections import OrderedDict, namedtuple
from collections.abc import Sequence
import itertools
import logging

import numpy as np
import scipy.linalg as la
import pandas as pd

from NNPDF import ThPredictions, CommonData, Experiment
from reportengine.checks import require_one, remove_outer, check_not_empty
from reportengine.table import table
from reportengine import collect

from validphys.checks import (check_cuts_considered, check_pdf_is_montecarlo,
                              check_speclabels_different, check_two_dataspecs,
                              check_dataset_cuts_match_theorycovmat,
                              check_experiment_cuts_match_theorycovmat)
from validphys.core import DataSetSpec, PDF, ExperimentSpec
from validphys.calcutils import (
    all_chi2, central_chi2, calc_chi2, calc_phi, bootstrap_values,
    get_df_block, regularize_covmat)

log = logging.getLogger(__name__)



class Result: pass


#TODO: Eventually,only one of (NNPDFDataResult, StatsResult) should survive
class NNPDFDataResult(Result):
    """A result fills its values from a libnnpf data object"""
    def __init__(self, dataobj):
        self._central_value = dataobj.get_cv()

    @property
    def central_value(self):
        return self._central_value

    def __len__(self):
        return len(self.central_value)

class StatsResult(Result):
    def __init__(self, stats):
        self.stats = stats

    @property
    def central_value(self):
        return self.stats.central_value()

    @property
    def std_error(self):
        return self.stats.std_error()

class DataResult(NNPDFDataResult):

    def __init__(self, dataobj, covmat, sqrtcovmat):
        super().__init__(dataobj)
        self._covmat = covmat
        self._sqrtcovmat = sqrtcovmat


    @property
    def label(self):
        return "Data"

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


class ThPredictionsResult(NNPDFDataResult):

    def __init__(self, dataobj, stats_class, label=None):
        self.stats_class = stats_class
        self.label = label
        self._std_error = dataobj.get_error()
        self._rawdata = dataobj.get_data()
        super().__init__(dataobj)

    @property
    def std_error(self):
        return self._std_error

    @staticmethod
    def make_label(pdf, dataset):
        """Deduce a reasonsble label for the result based on pdf and dataspec"""
        th = dataset.thspec
        if hasattr(pdf,'label'):
            if hasattr(th, 'label'):
                label = ' '.join((pdf.label, th.label))
            else:
                label = pdf.label
        elif hasattr(th, 'label'):
            label = th.label
        else:
            label = ('%s@<Theory %s>' % (pdf, th.id))
        return label


    @classmethod
    def from_convolution(cls, pdf, dataset, loaded_pdf=None, loaded_data=None):
        if loaded_pdf  is None:
            loaded_pdf = pdf.load()
        if loaded_data is None:
            loaded_data = dataset.load()
        th_predictions = ThPredictions(loaded_pdf, loaded_data)


        label = cls.make_label(pdf, dataset)


        return cls(th_predictions, pdf.stats_class, label)

class PositivityResult(StatsResult):
    @classmethod
    def from_convolution(cls, pdf, posset):
        loaded_pdf = pdf.load()
        loaded_pos = posset.load()
        data = loaded_pos.GetPredictions(loaded_pdf)
        stats = pdf.stats_class(data.T)
        return cls(stats)

    @property
    def rawdata(self):
        return self.stats.data

def experiments_index(experiments):
    """Return a pandas.MultiIndex with levels for
       experiment, dataset and point respectively"""
    records = []
    for experiment in experiments:
        for dataset in experiment.datasets:
            if dataset.cuts:
                ndata = len(dataset.cuts.load())
            else:
                #No cuts - use all data
                ndata = dataset.commondata.ndata
            for idat in range(ndata):
                records.append(
                    dict(
                        [('experiment', str(experiment.name)),
                         ('dataset', str(dataset.name)),
                         ('id', idat),]))

    columns = ['experiment', 'dataset', 'id']
    df = pd.DataFrame(records, columns=columns)
    df.set_index(columns, inplace=True)
    return df.index

def experiments_data(experiment_result_table):
    """Returns list of data values for the input experiments."""
    data_central_values = experiment_result_table["data_central"]
    return data_central_values

def experiment_result_table_no_table(experiments_results, experiments_index):
    """Generate a table containing the data central value, the central prediction,
    and the prediction for each PDF member."""
    result_records = []
    for experiment_results in experiments_results:
        dt, th = experiment_results
        for index, (dt_central, th_central) in enumerate(zip(dt.central_value, th.central_value)):
            replicas = (('rep_%05d'%(i+1), th_rep) for
                        i, th_rep in enumerate(th._rawdata[index, :]))

            result_records.append(dict([
                                 ('data_central', dt_central),
                                 ('theory_central', th_central),
                                  *replicas
                                 ]))
    if not result_records:
        log.warning("Empty records for experiment results")
        return pd.DataFrame()
    df =  pd.DataFrame(result_records, columns=result_records[0].keys(),
                       index=experiments_index)

    return df

@table
def experiment_result_table(experiment_result_table_no_table):
    """Duplicate of experiment_result_table_no_table but with a table decorator."""
    return experiment_result_table_no_table

@table
def experiment_result_table_68cl(experiment_result_table_no_table: pd.DataFrame, pdf: PDF):
    """Generate a table containing the data central value, the central prediction,
    and 68% confidence level bounds of the prediction.
    """
    df = experiment_result_table_no_table
    # replica data is every columns after central values, transpose for stats class
    replica_data = df.iloc[:, 2:].values.T
    # Use pdf stats class but reshape output to have each row as a data point
    stats = [level.reshape(-1, 1) for level in pdf.stats_class(replica_data).errorbar68()]
    # concatenate for dataframe construction
    stats_array = np.concatenate(stats, axis=1)
    df_cl = pd.DataFrame(
        stats_array,
        index=df.index,
        columns=['theory_lower', 'theory_upper'])
    res = pd.concat([df.iloc[:, :2], df_cl], axis=1)
    return res

experiments_covariance_matrix = collect('experiment_covariance_matrix', ('experiments',))


def experiments_covmat_no_table(
        experiments, experiments_index, experiments_covariance_matrix):
    """Export the covariance matrix for the experiments. It exports the full
    (symmetric) matrix, with the 3 first rows and columns being:

        - experiment name

        - dataset name

        - index of the point within the dataset.
    """
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, experiment_covmat in zip(
            experiments, experiments_covariance_matrix):
        name = experiment.name
        df.loc[[name],[name]] = experiment_covmat
    return df

@table
def experiments_covmat(experiments_covmat_no_table):
    """Duplicate of experiments_covmat_no_table but with a table decorator."""
    return experiments_covmat_no_table

exps_sqrt_covmat = collect(
    'experiment_sqrt_covariance_matrix',
    ('experiments',)
)

@table
def experiments_sqrtcovmat(
        experiments, experiments_index, exps_sqrt_covmat):
    """Like experiments_covmat, but dump the lower triangular part of the
    Cholesky decomposition as used in the fit. The upper part indices are set
    to zero.
    """
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, exp_sqrt_covmat in zip(
            experiments, exps_sqrt_covmat):
        name = experiment.name
        exp_sqrt_covmat[np.triu_indices_from(exp_sqrt_covmat, k=1)] = 0
        df.loc[[name],[name]] = exp_sqrt_covmat
    return df

@table
def experiments_invcovmat(
        experiments, experiments_index, experiments_covariance_matrix):
    """Compute and export the inverse covariance matrix.
    Note that this inverts the matrices with the LU method which is
    suboptimal."""
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment, experiment_covmat in zip(
            experiments, experiments_covariance_matrix):
        name = experiment.name
        #Improve this inversion if this method tuns out to be important
        invcov = la.inv(experiment_covmat)
        df.loc[[name],[name]] = invcov
    return df


@table
def experiments_normcovmat(experiments_covmat, experiments_data):
    """Calculates the experimental covariance matrix normalised to data."""
    df = experiments_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experiments_corrmat(experiments_covmat):
    """Generates the experimental correlation matrix with experiments_covmat as input"""
    df = experiments_covmat
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
    return mat

@table
def closure_pseudodata_replicas(experiments, pdf, nclosure:int,
                                experiments_index, nnoisy:int=0):
    """Generate closure pseudodata replicas from the given pdf.

    nclosure: Number of Level 1 pseudodata replicas.

    nnoisy:   Number of Level 2 replicas generated out of each pseudodata replica.

    The columns of the table are of the form (clos_0, noise_0_n0 ..., clos_1, ...)
    """

    #TODO: Do this somewhere else
    from  NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0,0)
    data = np.zeros((len(experiments_index), nclosure*(1+nnoisy)))

    cols = []
    for i in range(nclosure):
        cols += ['clos_%04d'%i, *['noise_%04d_%04d'%(i,j) for j in range(nnoisy)]]


    loaded_pdf = pdf.load()

    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        closure_exp = Experiment(exp.load())

        #TODO: This is probably computed somewhere else... All this code is
        #very error prone.
        #The predictions are for the unmodified experiment.
        predictions = [ThPredictions(loaded_pdf, d.load()) for d in exp]


        exp_location = experiments_index.get_loc(closure_exp.GetExpName())

        index = itertools.count()
        for i in range(nclosure):
            #Generate predictions with experimental noise, a different for
            #each closure set.
            closure_exp.MakeClosure(predictions, True)
            data[exp_location, next(index)] = closure_exp.get_cv()
            for j in range(nnoisy):
                #If we don't copy, we generate noise on top of the noise,
                #which is not what we want.
                replica_exp = Experiment(closure_exp)
                replica_exp.MakeReplica()

                data[exp_location, next(index)] = replica_exp.get_cv()


    df = pd.DataFrame(data, index=experiments_index,
                      columns=cols)

    return df


@check_dataset_cuts_match_theorycovmat
def covariance_matrix(
    dataset:DataSetSpec,
    fitthcovmat,
    t0set:(PDF, type(None)) = None,
    perform_covmat_reg=False,
    norm_threshold=3):
    """Returns the covariance matrix for a given `dataset`. By default the
    data central values will be used to calculate the multiplicative contributions
    to the covariance matrix.

    The matrix can instead be constructed with
    the t0 proceedure if the user sets `use_t0` to True and gives a
    `t0pdfset`. In this case the central predictions from the `t0pdfset` will be
    used to calculate the multiplicative contributions to the covariance matrix.
    More information on the t0 procedure can be found here:
    https://arxiv.org/abs/0912.2276

    The user can specify `use_fit_thcovmat_if_present` to be True
    and provide a corresponding `fit`. If the theory covmat was used in the
    corresponding `fit` and the specfied `dataset` was used in the fit then
    the theory covariance matrix for this `dataset` will be added in quadrature
    to the experimental covariance matrix.

    Covariance matrix can be regularized according to
    `calcutils.regularize_covmat` if the user specifies `perform_covmat_reg` to
    be true. This algorithm sets a minimum threshold for eigenvalues that the
    corresponding correlation matrix can have to be:

    1/(norm_threshold)^2

    which has the effect of limiting the L2 norm of the inverse of the correlation
    matrix. By default norm_threshold is set to 3.

    Parameters
    ----------
    dataset : DataSetSpec
        object parsed from the `dataset_input` runcard key
    fitthcovmat: None or ThCovMatSpec
        None if either `use_thcovmat_if_present` is False or if no theory
        covariance matrix was used in the corresponding fit
    t0set: None or PDF
        None if `use_t0` is False or a PDF parsed from `t0pdfset` runcard key
    perform_covmat_reg: bool
        whether or not to regularize the covariance matrix
    norm_threshold: number
        threshold used to regularize covariance matrix

    Returns
    -------
    covariance_matrix : array
        a covariance matrix as a numpy array

    Examples
    --------

    >>> from validphys.api import API
    >>> inp = {
            'dataset_input': {'ATLASTTBARTOT'},
            'theoryid': 52,
            'use_cuts': 'no_cuts'
        }
    >>> cov = API.covariance_matrix(**inp)
    TODO: complete example
    """
    loaded_data = dataset.load()

    if t0set:
        #Copy data to avoid chaos
        loaded_data = type(loaded_data)(loaded_data)
        log.debug("Setting T0 predictions for %s" % dataset)
        loaded_data.SetT0(t0set.load_t0())

    covmat = loaded_data.get_covmat()
    if fitthcovmat:
        loaded_thcov = fitthcovmat.load()
        covmat += get_df_block(loaded_thcov, dataset.name, level=1)
    if perform_covmat_reg:
        covmat = regularize_covmat(
            covmat,
            norm_threshold=norm_threshold
        )
    return covmat

datasets_covariance_matrix = collect(
    'covariance_matrix',
    ('experiments', 'experiment',)
)

def sqrt_covariance_matrix(covariance_matrix: np.array):
    """Returns the lower-triangular Cholesky factor `covariance_matrix`

    Parameters
    ----------
    covariance_matrix: array
        a positive definite covariance matrix

    Returns
    -------
    sqrt_covariance_matrix : array
        lower triangular Cholesky factor of covariance_matrix

    Notes
    -----
    The lower triangular is useful for efficient calculation of the chi^2

    Examples
    --------

    >>> import numpy as np
    >>> from validphys.results import sqrt_covariance_matrix
    >>> a = np.array([[1, 0.5], [0.5, 1]])
    >>> sqrt_covariance_matrix(a)
    array([[1.  , 0.    ],
        [0.5    , 0.8660254]])
    """
    return la.cholesky(covariance_matrix, lower=True)

@check_experiment_cuts_match_theorycovmat
def experiment_covariance_matrix(
        experiment: ExperimentSpec,
        fitthcovmat,
        t0set:(PDF, type(None)) = None,
        perform_covmat_reg=False,
        norm_threshold=3):
    """Like `covariance_matrix` except for an experiment"""
    loaded_data = experiment.load()

    if t0set:
        #Copy data to avoid chaos
        loaded_data = type(loaded_data)(loaded_data)
        log.debug("Setting T0 predictions for %s" % experiment)
        loaded_data.SetT0(t0set.load_t0())

    covmat = loaded_data.get_covmat()

    if fitthcovmat:
        loaded_thcov = fitthcovmat.load()
        ds_names = loaded_thcov.index.get_level_values(1)
        indices = np.in1d(ds_names, [ds.name for ds in experiment.datasets]).nonzero()[0]
        covmat += loaded_thcov.iloc[indices, indices].values
    if perform_covmat_reg:
        covmat = regularize_covmat(
            covmat,
            norm_threshold=norm_threshold
        )
    return covmat

def experiment_sqrt_covariance_matrix(experiment_covariance_matrix):
    """Like `sqrt_covariance_matrix` but for an experiment"""
    return sqrt_covariance_matrix(experiment_covariance_matrix)

def results(
        dataset:(DataSetSpec),
        pdf:PDF,
        covariance_matrix,
        sqrt_covariance_matrix
    ):
    """Tuple of data and theory results for a single pdf. The data will have an associated
    covariance matrix, which can include a contribution from the theory covariance matrix which
    is constructed from scale variation. The inclusion of this covariance matrix by default is used
    where available, however this behaviour can be modified with the flag `use_theorycovmat`.

    The theory is specified as part of the dataset.
    An experiment is also allowed.
    (as a result of the C++ code layout)."""
    data = dataset.load()
    return (DataResult(data, covariance_matrix, sqrt_covariance_matrix),
            ThPredictionsResult.from_convolution(pdf, dataset, loaded_data=data))

def experiment_results(
        experiment,
        pdf:PDF,
        experiment_covariance_matrix,
        experiment_sqrt_covariance_matrix):
    """Like `results` but for a whole experiment"""
    return results(
        experiment,
        pdf,
        experiment_covariance_matrix,
        experiment_sqrt_covariance_matrix
        )

#It's better to duplicate a few lines than to complicate the logic of
#``results`` to support this.
#TODO: The above comment doesn't make sense after adding T0. Deprecate this
def pdf_results(
        dataset:(DataSetSpec,  ExperimentSpec),
        pdfs:Sequence,
        covariance_matrix,
        sqrt_covariance_matrix):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    data = dataset.load()
    th_results = []
    for pdf in pdfs:
        th_result = ThPredictionsResult.from_convolution(pdf, dataset,
                                                         loaded_data=data)
        th_results.append(th_result)


    return (DataResult(data, covariance_matrix, sqrt_covariance_matrix), *th_results)

@require_one('pdfs', 'pdf')
@remove_outer('pdfs', 'pdf')
def one_or_more_results(dataset:(DataSetSpec, ExperimentSpec),
                        covariance_matrix,
                        sqrt_covariance_matrix,
                        pdfs:(type(None), Sequence)=None,
                        pdf:(type(None), PDF)=None):
    """Generate a list of results, where the first element is the data values,
    and the next is either the prediction for pdf or for each of the pdfs.
    Which of the two is selected intelligently depending on the namespace,
    when executing as an action."""
    if pdf:
        return results(dataset, pdf, covariance_matrix, sqrt_covariance_matrix)
    else:
        return pdf_results(dataset, pdfs, covariance_matrix, sqrt_covariance_matrix)
    raise ValueError("Either 'pdf' or 'pdfs' is required")


Chi2Data = namedtuple('Chi2Data', ('replica_result', 'central_result', 'ndata'))

def abs_chi2_data(results):
    """Return a tuple (member_chi², central_chi², numpoints) for a
    given dataset"""
    data_result, th_result = results

    chi2s = all_chi2(results)

    central_result = central_chi2(results)

    return Chi2Data(th_result.stats_class(chi2s[:, np.newaxis]),
                    central_result, len(data_result))


def abs_chi2_data_experiment(experiment_results):
    """Like `abs_chi2_data` but for a whole experiment"""
    return abs_chi2_data(experiment_results)

def phi_data(abs_chi2_data):
    """Calculate phi using values returned by `abs_chi2_data`.

    Returns tuple of (phi, numpoints)

    For more information on how phi is calculated see Eq.(24) in
    1410.8849
    """
    alldata, central, npoints = abs_chi2_data
    return (np.sqrt((alldata.data.mean() - central)/npoints), npoints)

def phi_data_experiment(abs_chi2_data_experiment):
    """Like `phi_data` but for whole experiment"""
    return phi_data(abs_chi2_data_experiment)

@check_pdf_is_montecarlo
def bootstrap_phi_data_experiment(experiment_results, bootstrap_samples=500):
    """Takes the data result and theory prediction for a given experiment and
    then returns a bootstrap distribution of phi.
    By default `bootstrap_samples` is set to a sensible value (500). However
    a different value can be specified in the runcard.

    For more information on how phi is calculated see `phi_data`
    """
    dt, th = experiment_results
    diff = np.array(th._rawdata - dt.central_value[:, np.newaxis])
    phi_resample = bootstrap_values(diff, bootstrap_samples,
                                    apply_func=(lambda x, y: calc_phi(y, x)),
                                    args=[dt.sqrtcovmat])
    return phi_resample

@check_pdf_is_montecarlo
def bootstrap_chi2_central_experiment(experiment_results, bootstrap_samples=500,
                                      boot_seed=123):
    """Takes the data result and theory prediction for a given experiment and
    then returns a bootstrap distribution of central chi2.
    By default `bootstrap_samples` is set to a sensible value (500). However
    a different value can be specified in the runcard.
    """
    dt, th = experiment_results
    diff = np.array(th._rawdata - dt.central_value[:, np.newaxis])
    cchi2 = lambda x, y: calc_chi2(y, x.mean(axis=1))
    chi2_central_resample = bootstrap_values(diff, bootstrap_samples, boot_seed=boot_seed,
                                             apply_func=(cchi2), args=[dt.sqrtcovmat])
    return chi2_central_resample

def chi2_breakdown_by_dataset(experiment_results, experiment, t0set,
                              prepend_total:bool=True,
                              datasets_sqrtcovmat=None) -> dict:
    """Return a dict with the central chi² of each dataset in the experiment,
    by breaking down the experiment results. If ``prepend_total`` is True.
    """
    dt, th = experiment_results
    sqrtcovmat = dt.sqrtcovmat
    central_diff = th.central_value - dt.central_value
    d = {}
    if prepend_total:
        d['Total'] = (calc_chi2(sqrtcovmat, central_diff), len(sqrtcovmat))


    #Allow lower level access useful for pseudodata and such.
    #TODO: This is a hack and we should get rid of it.
    if isinstance(experiment, Experiment):
        loaded_exp = experiment
    else:
        loaded_exp = experiment.load()

    #TODO: This is horrible. find a better way to do it.
    if t0set:
        loaded_exp = type(loaded_exp)(loaded_exp)
        loaded_exp.SetT0(t0set.load_T0())

    indmin = indmax = 0

    if datasets_sqrtcovmat is None:
        datasets_sqrtcovmat = (ds.get_sqrtcovmat() for ds in loaded_exp.DataSets())

    for ds, mat  in zip(loaded_exp.DataSets(), datasets_sqrtcovmat):
        indmax += len(ds)
        d[ds.GetSetName()] = (calc_chi2(mat, central_diff[indmin:indmax]), len(mat))
        indmin = indmax
    return d

def _chs_per_replica(chs):
    th, _, l = chs
    return th.data.ravel()/l


@table
def experiments_chi2_table(experiments, pdf, experiments_chi2,
                           each_dataset_chi2):
    """Return a table with the chi² to the experiments and each dataset on
    the experiments."""
    dschi2 = iter(each_dataset_chi2)
    records = []
    for experiment, expres in zip(experiments, experiments_chi2):
        stats = chi2_stats(expres)
        stats['experiment'] = experiment.name
        records.append(stats)
        for dataset, dsres in zip(experiment, dschi2):
            stats = chi2_stats(dsres)
            stats['experiment'] = dataset.name
            records.append(stats)
    return pd.DataFrame(records)

@check_cuts_considered
@table
def closure_shifts(experiments_index, fit, use_cuts, experiments):
    """Save the differenve between the fitted data and the real commondata
    values.

    Actually shifts is what should be saved in the first place, rather than
    thi confusing fiddling with Commondata, but until we can implement this at
    the C++ level, we just dave it here.
    """
    name, fitpath = fit
    result = np.zeros(len(experiments_index))
    for experiment in experiments:
        for dataset in experiment:
            dspath = fitpath/'filter'/dataset.name
            cdpath = dspath/("DATA_" + dataset.name + ".dat")
            try:
                syspath = next( (dspath/'systypes').glob('*.dat'))
            except StopIteration as e:
                raise FileNotFoundError("No systype "
                "file found in filter folder %s" % (dspath/'systypes')) from e
            cd = CommonData.ReadFile(str(cdpath), str(syspath))
            loc = experiments_index.get_loc((experiment.name, dataset.name))
            result[loc] = cd.get_cv() - dataset.load().get_cv()
    return pd.DataFrame(result, index=experiments_index)




def positivity_predictions(pdf, posdataset):
    """Return an object containing the values of the positivuty observable."""
    return PositivityResult.from_convolution(pdf, posdataset)

positivity_predictions_for_pdfs = collect(positivity_predictions, ('pdfs',))

def count_negative_points(possets_predictions):
    """Return the number of replicas with negative predictions for each bin
    in the positivity observable."""
    return np.sum([(r.rawdata < 0).sum(axis=1) for r in possets_predictions], axis=0)


chi2_stat_labels = {
    'central_mean': r'$<\chi^2_{0}>_{data}$',
    'npoints': r'$N_{data}$',
    'perreplica_mean': r'$\left< \chi^2 \right>_{rep,data}$',
    'perreplica_std': r'$\left<std_{rep}(\chi^2)\right>_{data}$',
    'chi2_per_data': r'$\frac{\chi^2}{N_{data}}$'
}

def chi2_stats(abs_chi2_data):
    """Compute severa estimators from the chi²:

     - central_mean

     - npoints

     - perreplica_mean

     - perreplica_std

     - chi2_per_data
    """
    rep_data, central_result, npoints = abs_chi2_data
    m = central_result.mean()
    rep_mean = rep_data.central_value().mean()
    return OrderedDict([
            ('central_mean'        ,  m),
            ('npoints'             , npoints),
            ('chi2_per_data', m/npoints),
            ('perreplica_mean', rep_mean),
            ('perreplica_std',  rep_data.std_error().mean()),
           ])


@table
def dataset_chi2_table(chi2_stats, dataset):
    """Show the chi² estimators for a given dataset"""
    return pd.DataFrame(chi2_stats, index=[dataset.name])

fits_experiment_chi2_data = collect(
    'experiments_chi2', ('fits', 'fit_context_groupby_experiment',))
fits_experiments = collect(
    'experiments', ('fits', 'fit_context_groupby_experiment',))

def fit_name_with_covmat_label(fit, fitthcovmat):
    """If theory covariance matrix is being used to calculate statistical estimators for the `fit`
    then appends (exp + th) onto the fit name for use in legends and column headers to help the user
    see what covariance matrix was used to produce the plot or table they are looking at.
    """
    if fitthcovmat:
        label = str(fit) + " (exp + th)"
    else:
        label = str(fit)
    return label

fits_name_with_covmat_label = collect('fit_name_with_covmat_label', ('fits',))


#TODO: Possibly get rid of the per_point_data parameter and have separate
#actions for absolute and relative tables.
@table
def fits_experiments_chi2_table(fits_name_with_covmat_label, fits_experiments,
                                fits_experiment_chi2_data, per_point_data:bool=True):
    """A table with the chi2 for each included experiment in the fits,
    computed with the theory corresponding to each fit.  If points_per_data is
    True, the chi² will be shown divided by ndata.
    Otherwise they will be absolute."""
    dfs = []
    cols = ('ndata', r'$\chi^2/ndata$') if per_point_data else ('ndata', r'$\chi^2$')
    for label, experiments, exps_chi2 in zip(
            fits_name_with_covmat_label, fits_experiments, fits_experiment_chi2_data):
        records = []
        for experiment, exp_chi2 in zip(experiments, exps_chi2):
            mean_chi2 = exp_chi2.central_result.mean()
            npoints = exp_chi2.ndata
            records.append(dict(
                experiment=str(experiment),
                npoints=npoints,
                mean_chi2 = mean_chi2

            ))
        df = pd.DataFrame.from_records(records,
                 columns=('experiment', 'npoints', 'mean_chi2'),
                 index = ('experiment', )
             )
        if per_point_data:
            df['mean_chi2'] /= df['npoints']
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    res =  pd.concat(dfs, axis=1)
    return res

fits_experiments_phi = collect(
    'experiments_phi', ('fits', 'fit_context_groupby_experiment'))

@table
def fits_experiments_phi_table(fits_name_with_covmat_label, fits_experiments, fits_experiments_phi):
    """For every fit, returns phi and number of data points per experiment, where experiment is
    a collection of datasets grouped according to the experiment key in the PLOTTING info file
    """
    dfs = []
    cols = ('ndata', r'$\phi$')
    for label, experiments, exps_phi in zip(
            fits_name_with_covmat_label, fits_experiments, fits_experiments_phi):
        records = []
        for experiment, (exp_phi, npoints) in zip(experiments, exps_phi):
            npoints = npoints
            records.append(dict(
                experiment=str(experiment),
                npoints=npoints,
                phi = exp_phi

            ))
        df = pd.DataFrame.from_records(records,
                 columns=('experiment', 'npoints', 'phi'),
                 index = ('experiment', )
             )
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    res =  pd.concat(dfs, axis=1)
    return res

@table
@check_speclabels_different
def dataspecs_experiments_chi2_table(dataspecs_speclabel, dataspecs_experiments,
                                     dataspecs_experiment_chi2_data,
                                     per_point_data:bool=True):
    """Same as fits_experiments_chi2_table but for an arbitrary list of dataspecs."""
    return fits_experiments_chi2_table(dataspecs_speclabel,
                                       dataspecs_experiments,
                                       dataspecs_experiment_chi2_data,
                                       per_point_data=per_point_data)


@table
def fits_datasets_chi2_table(fits_name_with_covmat_label, fits_experiments, fits_chi2_data,
                             per_point_data:bool=True):
    """A table with the chi2 for each included dataset in the fits, computed
    with the theory corresponding to the fit. The result are indexed in two
    levels by experiment and dataset, where experiment is the grouping of datasets according to the
    `experiment` key in the PLOTTING info file.  If points_per_data is True, the chi² will be shown
    divided by ndata. Otherwise they will be absolute."""

    chi2_it = iter(fits_chi2_data)

    cols = ('ndata', r'$\chi^2/ndata$') if per_point_data else ('ndata', r'$\chi^2$')

    dfs = []
    for label, experiments in zip(fits_name_with_covmat_label, fits_experiments):
        records = []
        for experiment in experiments:
            for dataset, chi2 in zip(experiment.datasets, chi2_it):
                ndata = chi2.ndata

                records.append(dict(
                    experiment=str(experiment),
                    dataset=str(dataset),
                    npoints=ndata,
                    mean_chi2 = chi2.central_result.mean()
                ))


        df = pd.DataFrame.from_records(records,
                 columns=('experiment', 'dataset', 'npoints', 'mean_chi2'),
                 index = ('experiment', 'dataset')
             )
        if per_point_data:
            df['mean_chi2'] /= df['npoints']
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    return pd.concat(dfs, axis=1)

@table
@check_speclabels_different
def dataspecs_datasets_chi2_table(dataspecs_speclabel, dataspecs_experiments,
                                  dataspecs_chi2_data, per_point_data:bool=True):
    """Same as fits_datasets_chi2_table but for arbitrary dataspecs."""
    return fits_datasets_chi2_table(dataspecs_speclabel, dataspecs_experiments,
                                    dataspecs_chi2_data, per_point_data=per_point_data)


#NOTE: This will need to be changed to work with `data` when that gets added
fits_total_chi2_data = collect('total_experiments_chi2data', ('fits', 'fitcontext'))
dataspecs_total_chi2_data = collect('total_experiments_chi2data', ('dataspecs',))

#TODO: Decide what to do with the horrible totals code.
@table
def fits_chi2_table(
        fits_total_chi2_data,
        fits_datasets_chi2_table,
        fits_experiments_chi2_table,
        show_total:bool=False):
    """Show the chi² of each and number of points of each dataset and experiment
    of each fit, where experiment is a group of datasets according to the `experiment` key in
    the PLOTTING info file, computed with the theory corresponding to the fit. Dataset that are not
    included in some fit appear as `NaN`
    """
    lvs = fits_experiments_chi2_table.index
    #The explicit call to list is because pandas gets confused otherwise
    expanded_index = pd.MultiIndex.from_product((list(lvs), ["Total"]))
    edf = fits_experiments_chi2_table.set_index(expanded_index)
    ddf = fits_datasets_chi2_table
    dfs = []
    #TODO: Better way to do the merge preserving the order?
    for lv in lvs:
        dfs.append(pd.concat((edf.loc[lv],ddf.loc[lv]), copy=False, axis=0))
    if show_total:
        total_points = np.array(
            [total_chi2_data.ndata for total_chi2_data in fits_total_chi2_data])
        total_chi = np.array(
            [total_chi2_data.central_result for total_chi2_data in fits_total_chi2_data])
        total_chi /= total_points
        row = np.zeros(len(total_points)*2)
        row[::2] = total_points
        row[1::2] = total_chi
        df = pd.DataFrame(np.atleast_2d(row),
                          columns=fits_experiments_chi2_table.columns,
                          index=['Total'])
        dfs.append(df)
        keys = [*lvs, 'Total']
    else:
        keys = lvs

    res = pd.concat(dfs, axis=0, keys=keys)
    return res

@table
def dataspecs_chi2_table(
    dataspecs_total_chi2_data,
    dataspecs_datasets_chi2_table,
    dataspecs_experiments_chi2_table,
    show_total: bool = False,
):
    """Same as fits_chi2_table but for an arbitrary list of dataspecs"""
    return fits_chi2_table(
        dataspecs_total_chi2_data,
        dataspecs_datasets_chi2_table,
        dataspecs_experiments_chi2_table,
        show_total,
    )


@table
@check_two_dataspecs
def dataspecs_chi2_differences_table(dataspecs, dataspecs_chi2_table):
    """Given two dataspecs, print the chi² (using dataspecs_chi2_table)
    and the difference between the first and the second."""
    df = dataspecs_chi2_table.copy()
    #TODO: Make this mind the number of points somehow
    diff = df.iloc[:,1] - df.iloc[:,3]
    df['difference'] = diff
    return df


def total_experiments_chi2data(pdf: PDF, experiments_chi2):
    """Return the central and member chi² values for the sum
    over all experiments. Values are not normalised to ndat.
    """
    ndata = 0
    central_chi2 = 0
    nmembers = len(experiments_chi2[0].replica_result.error_members())  # ugh
    member_chi2  = np.zeros(nmembers)

    for cd in experiments_chi2:
        # not sure why the transpose or [0] are needed here
        member_chi2  += cd.replica_result.error_members().T[0]
        central_chi2 += cd.central_result
        ndata += cd.ndata
    return Chi2Data(pdf.stats_class(member_chi2), central_chi2, ndata)


def total_experiments_chi2(total_experiments_chi2data):
    """Return the total chi²/ndata for the combination of all
    experiments."""
    return total_experiments_chi2data.central_result/total_experiments_chi2data.ndata

@table
@check_not_empty('experiments')
def perreplica_chi2_table(experiments, experiments_chi2):
    """Chi² per point for each replica for each experiment.
    Also outputs the total chi² per replica.
    The columns come in two levels: The first is the name of the experiment,
    and the second is the number of points."""

    chs = experiments_chi2
    total_chis = np.zeros((len(experiments) + 1, 1+ len(chs[0].replica_result.error_members())))
    ls = []
    for i,ch in enumerate(chs, 1):
        th, central, l = ch
        total_chis[i]= [central, *th.error_members()]
        ls.append(l)
    #total_chis/=total_l
    total_chis[0] = np.sum(total_chis[1:,:], axis=0)
    total_n = np.sum(ls)
    total_chis[0]/= total_n
    total_chis[1:,:]/= np.array(ls)[:, np.newaxis]

    columns = pd.MultiIndex.from_arrays(
            (['Total', *[str(exp) for exp in experiments]],
             [total_n, *ls]), names=['name', 'npoints'])
    return pd.DataFrame(total_chis.T, columns =columns)

@table
def theory_description(theoryid):
    """A table with the theory settings."""
    return pd.DataFrame(pd.Series(theoryid.get_description()), columns=[theoryid])

def experiments_central_values_no_table(experiment_result_table_no_table):
    """Returns a theoryid-dependent list of central theory predictions
    for a given experiment."""
    central_theory_values = experiment_result_table_no_table["theory_central"]
    return central_theory_values

@table
def experiments_central_values(experiment_result_table):
    """Duplicate of experiments_central_values_no_table but takes
    experiment_result_table rather than experiments_central_values_no_table,
    and has a table decorator."""
    central_theory_values = experiment_result_table["theory_central"]
    return central_theory_values

dataspecs_each_dataset_chi2 = collect("each_dataset_chi2", ("dataspecs",))
each_dataset = collect("dataset", ("experiments", "experiment"))
dataspecs_each_dataset = collect("each_dataset", ("dataspecs",))

@table
@check_speclabels_different
def dataspecs_dataset_chi2_difference_table(
    dataspecs_each_dataset, dataspecs_each_dataset_chi2, dataspecs_speclabel):
    r"""Returns a table with difference between the chi2 and the expected chi2
    in units of the expected chi2 standard deviation, given by

        chi2_diff = (\chi2 - N)/sqrt(2N)

    for each dataset for each dataspec.

    """
    dfs = []
    cols = [r"$(\chi^2 - N)/\sqrt{2N}$"]
    for label, datasets, chi2s in zip(
        dataspecs_speclabel, dataspecs_each_dataset, dataspecs_each_dataset_chi2):
        records = []
        for dataset, chi2 in zip(datasets, chi2s):
            ndata = chi2.ndata

            records.append(dict(
                dataset=str(dataset),
                chi2_stat = (chi2.central_result.mean() - ndata)/np.sqrt(2*ndata)
            ))


        df = pd.DataFrame.from_records(records,
                columns=("dataset", "chi2_stat"),
                index = ("dataset",)
            )
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    return pd.concat(dfs, axis=1)

datasets_covmat_no_reg = collect(
    "covariance_matrix", ("experiments", "experiment", "no_covmat_reg"))
datasets_covmat_reg = collect(
    "covariance_matrix", ("experiments", "experiment", "do_covmat_reg"))

@table
def datasets_covariance_matrix_differences_table(
    each_dataset, datasets_covmat_no_reg, datasets_covmat_reg):
    """For each dataset calculate and tabulate two max differences upon
    regularization given a value for `norm_threshold`:

    - max relative difference to the diagonal of the covariance matrix (%)
    - max absolute difference to the correlation matrix of each covmat

    """
    records = []
    for ds, reg, noreg in zip(
        each_dataset, datasets_covmat_reg, datasets_covmat_no_reg):
        cov_diag_rel_diff = np.diag(reg)/np.diag(noreg)
        d_reg = np.sqrt(np.diag(reg))
        d_noreg = np.sqrt(np.diag(noreg))
        corr_reg = reg/d_reg[:, np.newaxis]/d_reg[np.newaxis, :]
        corr_noreg = noreg/d_noreg[:, np.newaxis]/d_noreg[np.newaxis, :]
        corr_abs_diff = abs(corr_reg - corr_noreg)
        records.append(dict(
                dataset=str(ds),
                covdiff= np.max(abs(cov_diag_rel_diff- 1))*100, #make percentage
                corrdiff=np.max(corr_abs_diff)
            ))
    df = pd.DataFrame.from_records(records,
        columns=("dataset", "covdiff", "corrdiff"),
        index = ("dataset",)
        )
    df.columns = ["Variance rel. diff. (%)", "Correlation max abs. diff."]
    return df

dataspecs_covmat_diff_tables = collect(
    "datasets_covariance_matrix_differences_table", ("dataspecs",))

@check_speclabels_different
@table
def dataspecs_datasets_covmat_differences_table(
    dataspecs_speclabel, dataspecs_covmat_diff_tables
):
    """For each dataspec calculate and tabulate the two covmat differences
    described in `datasets_covariance_matrix_differences_table`
    (max relative difference in variance and max absolute correlation difference)

    """
    df = pd.concat( dataspecs_covmat_diff_tables, axis=1)
    cols = df.columns.get_level_values(0).unique()
    df.columns = pd.MultiIndex.from_product((dataspecs_speclabel, cols))
    return df

experiments_chi2 = collect(abs_chi2_data_experiment, ('experiments',))
each_dataset_chi2 = collect(abs_chi2_data, ('experiments', 'experiment'))

experiments_results = collect(experiment_results, ('experiments',))

experiments_phi = collect(phi_data_experiment, ('experiments',))

pdfs_total_chi2 = collect(total_experiments_chi2, ('pdfs',))

experiments_bootstrap_phi = collect(bootstrap_phi_data_experiment, ('experiments',))

experiments_bootstrap_chi2_central = collect(bootstrap_chi2_central_experiment,
                                             ('experiments',))

#These are convenient ways to iterate and extract varios data from fits
fits_chi2_data = collect(abs_chi2_data, ('fits', 'fitcontext', 'experiments', 'experiment'))

fits_total_chi2 = collect('total_experiments_chi2', ('fits', 'fitcontext'))

fits_total_chi2_for_experiments = collect('total_experiment_chi2',
                                          ('fits', 'fittheoryandpdf',
                                           'expspec', 'experiment'))

fits_pdf = collect('pdf', ('fits', 'fitpdf'))

#Dataspec is so
dataspecs_results = collect('results', ('dataspecs',))
dataspecs_chi2_data = collect(abs_chi2_data, ('dataspecs', 'experiments', 'experiment'))
dataspecs_experiment_chi2_data = collect('experiments_chi2', ('dataspecs',))
dataspecs_total_chi2 = collect('total_experiments_chi2', ('dataspecs',))
dataspecs_experiments_bootstrap_phi = collect('experiments_bootstrap_phi', ('dataspecs',))

dataspecs_speclabel = collect('speclabel', ('dataspecs',), element_default=None)
dataspecs_cuts = collect('cuts', ('dataspecs',))
dataspecs_experiments = collect('experiments', ('dataspecs',))
dataspecs_dataset = collect('dataset', ('dataspecs',))
dataspecs_commondata = collect('commondata', ('dataspecs',))
dataspecs_pdf = collect('pdf', ('dataspecs',))
dataspecs_fit = collect('fit', ('dataspecs',))
