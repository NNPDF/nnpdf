"""Module for handling logic and manipulation of covariance and correlation
matrices on different levels of abstraction
"""
import logging

import numpy as np
import pandas as pd
import scipy.linalg as la

from reportengine import collect
from reportengine.table import table
from validphys.calcutils import get_df_block, regularize_covmat
from validphys.checks import (
    check_cuts_considered,
    check_data_cuts_match_theorycovmat,
    check_dataset_cuts_match_theorycovmat,
    check_norm_threshold,
    check_pdf_is_montecarlo_or_hessian,
    check_speclabels_different,
)
from validphys.commondata import loaded_commondata_with_cuts
from validphys.convolution import central_predictions
from validphys.core import PDF, DataGroupSpec, DataSetSpec
from validphys.covmats_utils import construct_covmat, systematics_matrix
from validphys.results import ThPredictionsResult

log = logging.getLogger(__name__)

INTRA_DATASET_SYS_NAME = ("UNCORR", "CORR", "THEORYUNCORR", "THEORYCORR")


def covmat_from_systematics(
    loaded_commondata_with_cuts,
    dataset_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    _central_values=None,
):
    """Take the statistical uncertainty and systematics table from
    a :py:class:`validphys.coredata.CommonData` object and
    construct the covariance matrix accounting for correlations between
    systematics.

    If the systematic has the name ``SKIP`` then it is ignored in the
    construction of the covariance matrix.

    ADDitive or MULTiplicative systypes are handled by either multiplying
    the additive or multiplicative uncertainties respectively. We convert
    uncertainties so that they are all in the same units as the data:
        - Additive (ADD) systematics are left unchanged
        - multiplicative (MULT) systematics need to be converted from a
        percentage by multiplying by the central value
        and dividing by 100.

    Finally, the systematics are split into the five possible archetypes
    of systematic uncertainties: uncorrelated (UNCORR), correlated (CORR),
    theory uncorrelated (THEORYUNCORR), theory correlated (THEORYCORR) and
    special correlated (SPECIALCORR) systematics.

    Uncorrelated contributions from statistical error, uncorrelated and
    theory uncorrelated are added in quadrature to the diagonal of the covmat.

    The contribution to the covariance matrix arising due to
    correlated systematics is schematically ``A_correlated @ A_correlated.T``,
    where A_correlated is a matrix N_dat by N_sys. The total contribution
    from correlated systematics is found by adding together the result of
    mutiplying each correlated systematic matrix by its transpose
    (correlated, theory_correlated and special_correlated).

    For more information on the generation of the covariance matrix see the
    `paper <https://arxiv.org/pdf/hep-ph/0501067.pdf>`_
    outlining the procedure, specifically equation 2 and surrounding text.

    Parameters
    ----------

    loaded_commondata_with_cuts : validphys.coredata.CommonData
        CommonData which stores information about systematic errors,
        their treatment and description.
    dataset_input: validphys.core.DataSetInput
        Dataset settings, contains the weight for the current dataset.
        The returned covmat will be divided by the dataset weight if
        ``use_weights_in_covmat``. The default weight is 1, which means
        the returned covmat will be unmodified.
    use_weights_in_covmat: bool
        Whether to weight the covmat, True by default.
    norm_threshold: number
        threshold used to regularize covariance matrix
    _central_values : None, np.array
        1-D array containing alternative central values to combine with the
        multiplicative errors to calculate their absolute contributions. By
        default this is None, and the experimental central values are used. However, this
        can be used to calculate, for example, the t0 covariance matrix by
        using the predictions from the central member of the t0 pdf.

    Returns
    -------
    cov_mat: np.array
        Numpy array which is N_dat x N_dat (where N_dat is the number of data points after cuts)
        containing uncertainty and correlation information.

    Example
    -------
    In order to use this function, simply call it from the API

    >>> from validphys.api import API
    >>> inp = dict(
    ...     dataset_input={'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10},
    ...     theoryid=162,
    ...     use_cuts="internal"
    ... )
    >>> cov = API.covmat_from_systematics(**inp)
    >>> cov.shape
    (28, 28)

    """
    covmat = construct_covmat(
        loaded_commondata_with_cuts.stat_errors.to_numpy(),
        loaded_commondata_with_cuts.systematic_errors(_central_values),
    )
    if use_weights_in_covmat:
        covmat = covmat / dataset_input.weight
    if norm_threshold is not None:
        covmat = regularize_covmat(covmat, norm_threshold=norm_threshold)
    return covmat


def dataset_inputs_covmat_from_systematics(
    dataset_inputs_loaded_cd_with_cuts,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    _list_of_central_values=None,
    _only_additive=False,
):
    """Given a list containing :py:class:`validphys.coredata.CommonData` s,
    construct the full covariance matrix.

    This is similar to :py:meth:`covmat_from_systematics`
    except that special corr systematics are concatenated across all datasets
    before being multiplied by their transpose to give off block-diagonal
    contributions. The other systematics contribute to the block diagonal in the
    same way as :py:meth:`covmat_from_systematics`.

    Parameters
    ----------
    dataset_inputs_loaded_cd_with_cuts : list[validphys.coredata.CommonData]
        list of CommonData objects.
    data_input: list[validphys.core.DataSetInput]
        Settings for each dataset, each element contains the weight for the
        current dataset. The elements of the returned covmat for dataset
        i and j will be divided by sqrt(weight_i)*sqrt(weight_j), if
        ``use_weights_in_covmat``. The default weight is 1, which means
        the returned covmat will be unmodified.
    use_weights_in_covmat: bool
        Whether to weight the covmat, True by default.
    norm_threshold: number
        threshold used to regularize covariance matrix
    _list_of_central_values: None, list[np.array]
        list of 1-D arrays which contain alternative central values which are
        combined with the multiplicative errors to calculate their absolute
        contribution. By default this is None and the experimental central
        values are used.

    Returns
    -------
    cov_mat : np.array
        Numpy array which is N_dat x N_dat (where N_dat is the number of data points after cuts)
        containing uncertainty and correlation information.

    Example
    -------
    This function can be called directly from the API:

    >>> dsinps = [
    ...     {'dataset': 'NMC'},
    ...     {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
    ...     {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}
    ... ]
    >>> inp = dict(dataset_inputs=dsinps, theoryid=162, use_cuts="internal")
    >>> cov = API.dataset_inputs_covmat_from_systematics(**inp)
    >>> cov.shape
    (235, 235)

    Which properly accounts for all dataset settings and cuts.

    """
    special_corrs = []
    block_diags = []
    weights = []
    if _list_of_central_values is None:
        # want to just pass None to systematic_errors method
        _list_of_central_values = [None] * len(dataset_inputs_loaded_cd_with_cuts)

    for cd, dsinp, central_values in zip(
        dataset_inputs_loaded_cd_with_cuts, data_input, _list_of_central_values
    ):
        # used if we want to separate additive and multiplicative errors in make_replica
        if _only_additive:
            sys_errors = cd.additive_errors
        else:
            sys_errors = cd.systematic_errors(central_values)
        stat_errors = cd.stat_errors.to_numpy()
        weights.append(np.full_like(stat_errors, dsinp.weight))
        # separate out the special uncertainties which can be correlated across
        # datasets
        is_intra_dataset_error = sys_errors.columns.isin(INTRA_DATASET_SYS_NAME)
        block_diags.append(construct_covmat(stat_errors, sys_errors.loc[:, is_intra_dataset_error]))
        special_corrs.append(sys_errors.loc[:, ~is_intra_dataset_error])

    # concat systematics across datasets
    special_sys = pd.concat(special_corrs, axis=0, sort=False)
    # non-overlapping systematics are set to NaN by concat, fill with 0 instead.
    special_sys.fillna(0, inplace=True)

    diag = la.block_diag(*block_diags)
    covmat = diag + special_sys.to_numpy() @ special_sys.to_numpy().T
    if use_weights_in_covmat:
        # concatenate weights and sqrt
        sqrt_weights = np.sqrt(np.concatenate(weights))
        # returns C_ij / (sqrt(w_i) * sqrt(w_j))
        covmat = (covmat / sqrt_weights).T / sqrt_weights
    if norm_threshold is not None:
        covmat = regularize_covmat(covmat, norm_threshold=norm_threshold)
    return covmat


@check_cuts_considered
def dataset_t0_predictions(dataset, t0set):
    """Returns the t0 predictions for a ``dataset`` which are the predictions
    calculated using the central member of ``pdf``. Note that if ``pdf`` has
    errortype ``replicas``, and the dataset is a hadronic observable then the
    predictions of the central member are subtly different to the central
    value of the replica predictions.

    Parameters
    ----------
    dataset: validphys.core.DataSetSpec
        dataset for which to calculate t0 predictions
    t0set: validphys.core.PDF
        pdf used to calculate the predictions

    Returns
    -------
    t0_predictions: np.array
        1-D numpy array with predictions for each of the cut datapoints.

    """
    # reshape because the underlying data has shape ndata * 1
    # accounting for the fact that some datasets are single datapoint
    return central_predictions(dataset, t0set).to_numpy().reshape(-1)


def t0_covmat_from_systematics(
    loaded_commondata_with_cuts,
    *,
    dataset_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    dataset_t0_predictions,
):
    """Like :py:func:`covmat_from_systematics` except uses the t0 predictions
    to calculate the absolute constributions to the covmat from multiplicative
    uncertainties. For more info on the t0 predictions see
    :py:func:`validphys.commondata.dataset_t0_predictions`.

    Parameters
    ----------
    loaded_commondata_with_cuts: validphys.coredata.CommonData
        commondata object for which to generate the covmat.
    dataset_input: validphys.core.DataSetInput
        Dataset settings, contains the weight for the current dataset.
        The returned covmat will be divided by the dataset weight if
        ``use_weights_in_covmat``. The default weight is 1, which means
        the returned covmat will be unmodified.
    use_weights_in_covmat: bool
        Whether to weight the covmat, True by default.
    dataset_t0_predictions: np.array
        1-D array with t0 predictions.

    Returns
    -------
    t0_covmat: np.array
        t0 covariance matrix

    """
    return covmat_from_systematics(
        loaded_commondata_with_cuts,
        dataset_input,
        use_weights_in_covmat,
        norm_threshold=norm_threshold,
        _central_values=dataset_t0_predictions,
    )


dataset_inputs_t0_predictions = collect("dataset_t0_predictions", ("data",))


def dataset_inputs_t0_covmat_from_systematics(
    dataset_inputs_loaded_cd_with_cuts,
    *,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    dataset_inputs_t0_predictions,
):
    """Like :py:func:`t0_covmat_from_systematics` except for all data

    Parameters
    ----------
    dataset_inputs_loaded_cd_with_cuts: list[validphys.coredata.CommonData]
        The CommonData for all datasets defined in ``dataset_inputs``.
    data_input: list[validphys.core.DataSetInput]
        Settings for each dataset, each element contains the weight for the
        current dataset. The elements of the returned covmat for dataset
        i and j will be divided by sqrt(weight_i)*sqrt(weight_j), if
        ``use_weights_in_covmat``. The default weight is 1, which means
        the returned covmat will be unmodified.
    use_weights_in_covmat: bool
        Whether to weight the covmat, True by default.
    dataset_inputs_t0_predictions: list[np.array]
        The t0 predictions for all datasets.

    Returns
    -------
    t0_covmat: np.array
        t0 covariance matrix matrix for list of datasets.
    """
    return dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts,
        data_input,
        use_weights_in_covmat,
        norm_threshold=norm_threshold,
        _list_of_central_values=dataset_inputs_t0_predictions,
    )


def dataset_inputs_t0_total_covmat_separate(
    dataset_inputs_t0_exp_covmat_separate, loaded_theory_covmat
):
    """
    Function to compute the covmat to be used for the sampling by make_replica.
    In this case the t0 prescription is used for the experimental covmat and the multiplicative
    errors are separated. Moreover, the theory covmat is added to experimental covmat.
    """
    covmat = dataset_inputs_t0_exp_covmat_separate
    covmat += loaded_theory_covmat
    return covmat


def dataset_inputs_t0_exp_covmat_separate(
    dataset_inputs_loaded_cd_with_cuts,
    *,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    dataset_inputs_t0_predictions,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica.
    In this case the t0 prescription is used for the experimental covmat and the multiplicative
    errors are separated.
    """
    covmat = generate_exp_covmat(
        dataset_inputs_loaded_cd_with_cuts,
        data_input,
        use_weights_in_covmat,
        norm_threshold,
        dataset_inputs_t0_predictions,
        True,
    )
    return covmat


def dataset_inputs_total_covmat_separate(
    dataset_inputs_exp_covmat_separate,
    loaded_theory_covmat,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica.
    In this case the t0 prescription is not used for the experimental covmat and the multiplicative
    errors are separated. Moreover, the theory covmat is added to experimental covmat.
    """
    covmat = dataset_inputs_exp_covmat_separate
    covmat += loaded_theory_covmat
    return covmat


def dataset_inputs_exp_covmat_separate(
    dataset_inputs_loaded_cd_with_cuts,
    *,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica.
    In this case the t0 prescription is not used for the experimental covmat and the multiplicative
    errors are separated.
    """
    covmat = generate_exp_covmat(
        dataset_inputs_loaded_cd_with_cuts,
        data_input,
        use_weights_in_covmat,
        norm_threshold,
        None,
        True,
    )
    return covmat


def dataset_inputs_t0_total_covmat(
    dataset_inputs_t0_exp_covmat,
    loaded_theory_covmat,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica and for the chi2
    by fitting_data_dict. In this case the t0 prescription is used for the experimental covmat
    and the multiplicative errors are included in it. Moreover, the theory covmat is added to experimental covmat.
    """
    covmat = dataset_inputs_t0_exp_covmat
    covmat += loaded_theory_covmat
    return covmat


def dataset_inputs_t0_exp_covmat(
    dataset_inputs_loaded_cd_with_cuts,
    *,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
    dataset_inputs_t0_predictions,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica and for the chi2
    by fitting_data_dict. In this case the t0 prescription is used for the experimental covmat
    and the multiplicative errors are included in it.
    """
    covmat = generate_exp_covmat(
        dataset_inputs_loaded_cd_with_cuts,
        data_input,
        use_weights_in_covmat,
        norm_threshold,
        dataset_inputs_t0_predictions,
        False,
    )
    return covmat


def dataset_inputs_total_covmat(
    dataset_inputs_exp_covmat,
    loaded_theory_covmat,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica and for the chi2
    by fitting_data_dict. In this case the t0 prescription is not used for the experimental covmat
    and the multiplicative errors are included in it. Moreover, the theory covmat is added to experimental covmat.
    """
    covmat = dataset_inputs_exp_covmat
    covmat += loaded_theory_covmat
    return covmat


def dataset_inputs_exp_covmat(
    dataset_inputs_loaded_cd_with_cuts,
    *,
    data_input,
    use_weights_in_covmat=True,
    norm_threshold=None,
):
    """
    Function to compute the covmat to be used for the sampling by make_replica and for the chi2
    by fitting_data_dict. In this case the t0 prescription is not used for the experimental covmat
    and the multiplicative errors are included in it.
    """
    covmat = generate_exp_covmat(
        dataset_inputs_loaded_cd_with_cuts,
        data_input,
        use_weights_in_covmat,
        norm_threshold,
        None,
        False,
    )
    return covmat


def generate_exp_covmat(
    datasets_input, data, use_weights, norm_threshold, _list_of_c_values, only_add
):
    """
    Function to generate the experimental covmat eventually using the t0 prescription. It is also
    possible to compute it only with the additive errors.

    Parameters
    ----------
        dataset_inputs: list[validphys.coredata.CommonData]
            list of CommonData objects.
        data: list[validphys.core.DataSetInput]
            Settings for each dataset, each element contains the weight for the
            current dataset. The elements of the returned covmat for dataset
            i and j will be divided by sqrt(weight_i)*sqrt(weight_j), if
            ``use_weights_in_covmat``. The default weight is 1, which means
            the returned covmat will be unmodified.
        use_weights: bool
            Whether to weight the covmat, True by default.
        norm_threshold: number
            threshold used to regularize covariance matrix
        _list_of_c_values: None, list[np.array]
            list of 1-D arrays which contain alternative central values which are
            combined with the multiplicative errors to calculate their absolute
            contribution. By default this is None and the experimental central
            values are used.
        only_add: bool
            specifies whether to use only the additive errors to compute the covmat

    Returns
    -------
        : np.array
        experimental covariance matrix
    """
    return dataset_inputs_covmat_from_systematics(
        datasets_input,
        data,
        use_weights,
        norm_threshold=norm_threshold,
        _list_of_central_values=_list_of_c_values,
        _only_additive=only_add,
    )


def sqrt_covmat(covariance_matrix):
    """Function that computes the square root of the covariance matrix.

    Parameters
    ----------
    covariance_matrix : np.array
        A positive definite covariance matrix, which is N_dat x N_dat (where
        N_dat is the number of data points after cuts) containing uncertainty
        and correlation information.

    Returns
    -------
    sqrt_mat : np.array
        The square root of the input covariance matrix, which is N_dat x N_dat
        (where N_dat is the number of data points after cuts), and which is the
        the lower triangular decomposition. The following should be ``True``:
        ``np.allclose(sqrt_covmat @ sqrt_covmat.T, covariance_matrix)``.

    Notes
    -----
    The square root is found by using the Cholesky decomposition. However, rather
    than finding the decomposition of the covariance matrix directly, the (upper
    triangular) decomposition is found of the corresponding correlation matrix
    and then the output of this is rescaled and then transposed as
    ``sqrt_matrix = (decomp * sqrt_diags).T``, where ``decomp`` is the Cholesky
    decomposition of the correlation matrix and ``sqrt_diags`` is the square root
    of the diagonal entries of the covariance matrix. This method is useful in
    situations in which the covariance matrix is near-singular. See
    `here <https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition>`_
    for more discussion on this.

    The lower triangular is useful for efficient calculation of the :math:`\chi^2`

    Example
    -------
    >>> import numpy as np
    >>> from validphys.api import API
    >>> API.sqrt_covmat(dataset_input={"dataset":"NMC"}, theoryid=162, use_cuts="internal")
    array([[0.0326543 , 0.        , 0.        , ..., 0.        , 0.        ,
            0.        ],
        [0.00314523, 0.01467259, 0.        , ..., 0.        , 0.        ,
            0.        ],
        [0.0037817 , 0.00544256, 0.02874822, ..., 0.        , 0.        ,
            0.        ],
        ...,
        [0.00043404, 0.00031169, 0.00020489, ..., 0.00441073, 0.        ,
            0.        ],
        [0.00048717, 0.00033792, 0.00022971, ..., 0.00126704, 0.00435696,
            0.        ],
        [0.00067353, 0.00050372, 0.0003203 , ..., 0.00107255, 0.00065041,
            0.01002952]])
    >>> sqrt_cov = API.sqrt_covmat(dataset_input={"dataset":"NMC"}, theoryid=162, use_cuts="internal")
    >>> cov = API.covariance_matrix(dataset_input={"dataset":"NMC"}, theoryid=162, use_cuts="internal")
    >>> np.allclose(np.linalg.cholesky(cov), sqrt_cov)
    True

    """
    dimensions = covariance_matrix.shape

    if covariance_matrix.size == 0:
        return np.zeros((0,0))
    elif dimensions[0] != dimensions[1]:
        raise ValueError(
            "The input covariance matrix should be square but "
            f"instead it has dimensions {dimensions[0]} x "
            f"{dimensions[1]}"
        )

    sqrt_diags = np.sqrt(np.diag(covariance_matrix))
    correlation_matrix = covariance_matrix / sqrt_diags[:, np.newaxis] / sqrt_diags
    decomp = la.cholesky(correlation_matrix)
    sqrt_matrix = (decomp * sqrt_diags).T
    return sqrt_matrix


def groups_covmat_no_table(groups_data, groups_index, groups_covmat_collection):
    """Export the covariance matrix for the groups. It exports the full
    (symmetric) matrix, with the 3 first rows and columns being:

        - group name

        - dataset name

        - index of the point within the dataset.
    """
    data = np.zeros((len(groups_index), len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_covmat in zip(groups_data, groups_covmat_collection):
        name = group.name
        df.loc[[name], [name]] = group_covmat
    return df


@table
def groups_covmat(groups_covmat_no_table):
    """Duplicate of groups_covmat_no_table but with a table decorator."""
    return groups_covmat_no_table


@table
def groups_sqrtcovmat(groups_data, groups_index, groups_sqrt_covmat):
    """Like groups_covmat, but dump the lower triangular part of the
    Cholesky decomposition as used in the fit. The upper part indices are set
    to zero.
    """
    data = np.zeros((len(groups_index), len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_sqrt_covmat in zip(groups_data, groups_sqrt_covmat):
        name = group.name
        group_sqrt_covmat[np.triu_indices_from(group_sqrt_covmat, k=1)] = 0
        df.loc[[name], [name]] = group_sqrt_covmat
    return df


@table
def groups_invcovmat(groups_data, groups_index, groups_covmat_collection):
    """Compute and export the inverse covariance matrix.
    Note that this inverts the matrices with the LU method which is
    suboptimal."""
    data = np.zeros((len(groups_index), len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_covmat in zip(groups_data, groups_covmat_collection):
        name = group.name
        # Improve this inversion if this method tuns out to be important
        invcov = la.inv(group_covmat)
        df.loc[[name], [name]] = invcov
    return df


@table
def groups_normcovmat(groups_covmat, groups_data_values):
    """Calculates the grouped experimental covariance matrix normalised to data."""
    df = groups_covmat
    groups_data_array = np.array(groups_data_values)
    mat = df / np.outer(groups_data_array, groups_data_array)
    return mat


@table
def groups_corrmat(groups_covmat):
    """Generates the grouped experimental correlation matrix with groups_covmat as input"""
    df = groups_covmat
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat)) ** (-0.5)
    mat = diag_minus_half[:, np.newaxis] * df * diag_minus_half
    return mat


@check_pdf_is_montecarlo_or_hessian
def pdferr_plus_covmat(dataset, pdf, covmat_t0_considered):
    """For a given `dataset`, returns the sum of the covariance matrix given by
    `covmat_t0_considered` and the PDF error:
    - If the PDF error_type is 'replicas', a covariance matrix is estimated from
      the replica theory predictions
    - If the PDF error_type is 'symmhessian', a covariance matrix is estimated using
      formulas from (mc2hessian) https://arxiv.org/pdf/1505.06736.pdf
    - If the PDF error_type is 'hessian' a covariance matrix is estimated using
      the hessian formula from Eq. 5 of https://arxiv.org/pdf/1401.0013.pdf


    Parameters
    ----------
    dataset: DataSetSpec
        object parsed from the `dataset_input` runcard key
    pdf: PDF
        monte carlo pdf used to estimate PDF error
    covmat_t0_considered: np.array
        experimental covariance matrix with the t0 considered

    Returns
    -------
    covariance_matrix: np.array
        sum of the experimental and pdf error as a numpy array

    Examples
    --------

    `use_pdferr` makes this action be used for `covariance_matrix`

    >>> from validphys.api import API
    >>> from import numpy as np
    >>> inp = {
            'dataset_input': {'dataset' : 'ATLASTTBARTOT'},
            'theoryid': 53,
            'pdf': 'NNPDF31_nlo_as_0118',
            'use_cuts': 'nocuts'
        }
    >>> a = API.covariance_matrix(**inp, use_pdferr=True)
    >>> b = API.pdferr_plus_covmat(**inp)
    >>> np.allclose(a == b)
    True
    """
    th = ThPredictionsResult.from_convolution(pdf, dataset)

    if pdf.error_type == 'replicas':
        pdf_cov = np.cov(th.error_members, rowvar=True)

    elif pdf.error_type == 'symmhessian':
        rescale_fac = pdf._rescale_factor()
        hessian_eigenvectors = th.error_members
        central_predictions = th.central_value

        # need to subtract the central set which is not the same as the average of the
        # Hessian eigenvectors.
        X = hessian_eigenvectors - central_predictions.reshape((central_predictions.shape[0], 1))
        # need to rescale the Hessian eigenvectors in case the eigenvector confidence interval is not 68%
        X = X / rescale_fac
        pdf_cov = X @ X.T

    elif pdf.error_type == 'hessian':
        rescale_fac = pdf._rescale_factor()
        hessian_eigenvectors = th.error_members
        
        # see core.HessianStats
        X = (hessian_eigenvectors[:,0::2] - hessian_eigenvectors[:,1::2])*0.5
        # need to rescale the Hessian eigenvectors in case the eigenvector confidence interval is not 68%
        X = X / rescale_fac
        pdf_cov = X @ X.T

    return pdf_cov + covmat_t0_considered


def reorder_thcovmat_as_expcovmat(fitthcovmat, data):
    """
    Reorder the thcovmat in such a way to match the order of the experimental covmat, which
    means the order of the runcard
    """
    theory_covmat = fitthcovmat.load()
    bb = [str(i) for i in data]
    tmp = theory_covmat.droplevel(0, axis=0).droplevel(0, axis=1)
    return tmp.reindex(index=bb, columns=bb, level=0)


def pdferr_plus_dataset_inputs_covmat(data, pdf, dataset_inputs_covmat_t0_considered, fitthcovmat):
    """Like `pdferr_plus_covmat` except for an experiment"""
    # do checks get performed here?
    if fitthcovmat is not None:
        # change ordering according to exp_covmat (so according to runcard order)
        return pdferr_plus_covmat(
            data,
            pdf,
            dataset_inputs_covmat_t0_considered
            + reorder_thcovmat_as_expcovmat(fitthcovmat, data).values,
        )
    return pdferr_plus_covmat(data, pdf, dataset_inputs_covmat_t0_considered)


def dataset_inputs_sqrt_covmat(dataset_inputs_covariance_matrix):
    """Like `sqrt_covmat` but for an group of datasets"""
    return sqrt_covmat(dataset_inputs_covariance_matrix)


def systematics_matrix_from_commondata(
    loaded_commondata_with_cuts, dataset_input, use_weights_in_covmat=True, _central_values=None
):
    """Returns a systematics matrix, :math:`A`, for the corresponding dataset.
    The systematics matrix is a square root of the covmat:

    .. math::

        C = A A^T

    and is obtained by concatenating a block diagonal of the uncorrelated uncertainties
    with the correlated systematics.

    """
    sqrt_covmat = systematics_matrix(
        loaded_commondata_with_cuts.stat_errors.to_numpy(),
        loaded_commondata_with_cuts.systematic_errors(_central_values),
    )
    if use_weights_in_covmat:
        return sqrt_covmat / np.sqrt(dataset_input.weight)
    return sqrt_covmat


def covmat_stability_characteristic(systematics_matrix_from_commondata):
    """
    Return a number characterizing the stability of an experimental covariance
    matrix against uncertainties in the correlation. It is defined as the L2
    norm (largest singular value) of the square root of the inverse correlation
    matrix. This is equivalent to the square root of the inverse of the
    smallest singular value of the correlation matrix:

    Z = (1/λ⁰)^½

    Where λ⁰ is the smallest eigenvalue of the correlation matrix.

    This is the number used as
    threshold in :py:func:`calcutils.regularize_covmat`. The interpretation
    is roughly what precision does the worst correlation need to
    have in order to not affect meaningfully the χ² computed using the
    covariance matrix, so for example a stability characteristic of 4 means
    that correlations need to be known with uncetainties less than 0.25.

    Examples
    --------

    >>> from validphys.api import API
    >>> API.covmat_stability_characteristic(dataset_input={"dataset": "NMC"},
    ... theoryid=162, use_cuts="internal")
    2.742658604186114

    """
    sqrtcov = systematics_matrix_from_commondata
    # copied from calcutils.regularize_l2 but just return stability condition.
    d = np.sqrt(np.sum(sqrtcov**2, axis=1))[:, np.newaxis]
    sqrtcorr = sqrtcov / d
    _, s, _ = la.svd(sqrtcorr, full_matrices=False)
    return 1 / s[-1]


dataset_inputs_stability = collect('covmat_stability_characteristic', ('dataset_inputs',))


@table
def dataset_inputs_stability_table(dataset_inputs_stability, dataset_inputs):
    """Return a table with py:func:`covmat_stability_characteristic` for all
    dataset inputs"""
    res = {}
    for ds, stab in zip(dataset_inputs, dataset_inputs_stability):
        res[ds.name] = stab

    return pd.Series(res, name="stability").sort_values()


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


@table
@check_norm_threshold
def datasets_covmat_differences_table(
    each_dataset, datasets_covmat_no_reg, datasets_covmat_reg, norm_threshold
):
    """For each dataset calculate and tabulate two max differences upon
    regularization given a value for `norm_threshold`:

    - max relative difference to the diagonal of the covariance matrix (%)
    - max absolute difference to the correlation matrix of each covmat

    """
    records = []
    for ds, reg, noreg in zip(each_dataset, datasets_covmat_reg, datasets_covmat_no_reg):
        cov_diag_rel_diff = np.diag(reg) / np.diag(noreg)
        d_reg = np.sqrt(np.diag(reg))
        d_noreg = np.sqrt(np.diag(noreg))
        corr_reg = reg / d_reg[:, np.newaxis] / d_reg[np.newaxis, :]
        corr_noreg = noreg / d_noreg[:, np.newaxis] / d_noreg[np.newaxis, :]
        corr_abs_diff = abs(corr_reg - corr_noreg)
        records.append(
            dict(
                dataset=str(ds),
                covdiff=np.max(abs(cov_diag_rel_diff - 1)) * 100,  # make percentage
                corrdiff=np.max(corr_abs_diff),
            )
        )
    df = pd.DataFrame.from_records(
        records, columns=("dataset", "covdiff", "corrdiff"), index=("dataset",)
    )
    df.columns = ["Variance rel. diff. (%)", "Correlation max abs. diff."]
    return df


@check_speclabels_different
@table
def dataspecs_datasets_covmat_differences_table(dataspecs_speclabel, dataspecs_covmat_diff_tables):
    """For each dataspec calculate and tabulate the two covmat differences
    described in `datasets_covmat_differences_table`
    (max relative difference in variance and max absolute correlation difference)

    """
    df = pd.concat(dataspecs_covmat_diff_tables, axis=1)
    cols = df.columns.get_level_values(0).unique()
    df.columns = pd.MultiIndex.from_product((dataspecs_speclabel, cols))
    return df


def _covmat_t0_considered(covmat_t0_considered, fitthcovmat=None, dataset=None):
    """Helper function so we can dispatch the full
    covariance matrix, having considered both ``use_t0``
    and ``use_pdferr``
    """
    if fitthcovmat is not None:
        # exploit `reorder_thcovmat_as_expcovmat` to take only the part of the covmat for the relevant dataset
        return (
            covmat_t0_considered
            + reorder_thcovmat_as_expcovmat(fitthcovmat, [dataset]).values
        )
    return covmat_t0_considered

def _dataset_inputs_covmat_t0_considered(dataset_inputs_covmat_t0_considered, fitthcovmat, data):
    """Helper function so we can dispatch the full
    covariance matrix accross dataset_inputs, having considered both ``use_t0``
    and ``use_pdferr``
    """
    if fitthcovmat is not None:
        # change ordering according to exp_covmat (so according to runcard order)
        return (
            dataset_inputs_covmat_t0_considered
            + reorder_thcovmat_as_expcovmat(fitthcovmat, data).values
        )
    return dataset_inputs_covmat_t0_considered


groups_covmat_collection = collect(
    'dataset_inputs_covariance_matrix', ('group_dataset_inputs_by_metadata',)
)

groups_sqrt_covmat = collect('dataset_inputs_sqrt_covmat', ('group_dataset_inputs_by_metadata',))

dataspecs_covmat_diff_tables = collect("datasets_covmat_differences_table", ("dataspecs",))

fits_name_with_covmat_label = collect('fit_name_with_covmat_label', ('fits',))

datasets_covmat_no_reg = collect("covariance_matrix", ("data", "no_covmat_reg"))

datasets_covmat_reg = collect("covariance_matrix", ("data",))

datasets_covmat = collect('covariance_matrix', ('data',))

datasets_covariance_matrix = collect(
    'covariance_matrix',
    (
        'experiments',
        'experiment',
    ),
)
