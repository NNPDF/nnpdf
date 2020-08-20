"""Module for handling logic and manipulation of covariance and correlation
matrices on different levels of abstraction
"""
import logging

import numpy as np
import pandas as pd
import scipy.linalg as la

from reportengine import collect
from reportengine.table import table

from validphys.calcutils import regularize_covmat, get_df_block
from validphys.core import PDF, DataGroupSpec, DataSetSpec
from validphys.checks import (
    check_dataset_cuts_match_theorycovmat,
    check_norm_threshold,
    check_pdf_is_montecarlo,
    check_speclabels_different,
    check_data_cuts_match_theorycovmat,
)
from validphys.results import ThPredictionsResult

log = logging.getLogger(__name__)


def covmat_from_systematics(commondata):
    """Function that takes in a :py:class:`validphys.coredata.CommonData` object
    and generates the associated covariance matrix using the systematics. The logic
    is best described by the now deprecated C++ code:

    .. code-block:: c++

            auto CovMat = NNPDF::matrix<double>(ndat, ndat);

            for (int i = 0; i < ndat; i++)
            {
              for (int j = 0; j < ndat; j++)
              {
                double sig    = 0.0;
                double signor = 0.0;

                // Statistical error
                if (i == j)
                  sig += pow(stat_error[i],2);

                for (int l = 0; l < nsys; l++)
                {
                  sysError const& isys = systematic_errors[i][l];
                  sysError const& jsys = systematic_errors[j][l];
                  if (isys.name != jsys.name)
                      throw RuntimeException("ComputeCovMat", "Inconsistent naming of systematics");
                  if (isys.name == "SKIP")
                      continue; // Continue if systype is skipped
                  if ((isys.name == "THEORYCORR" || isys.name == "THEORYUNCORR") && !use_theory_errors)
                      continue; // Continue if systype is theoretical and use_theory_errors == false
                  const bool is_correlated = ( isys.name != "UNCORR" && isys.name !="THEORYUNCORR");
                  if (i == j || is_correlated)
                    switch (isys.type)
                    {
                      case ADD:   sig    += isys.add *jsys.add;  break;
                      case MULT: if (mult_errors) { signor += isys.mult*jsys.mult; break; }
                                 else { continue; }
                      case UNSET: throw RuntimeException("ComputeCovMat", "UNSET systype encountered");
                    }
                }

                // Covariance matrix entry
                CovMat(i, j) = (sig + signor*central_values[i]*central_values[j]*1e-4);
            // Covariance matrix weight
            CovMat(i, j) /= sqrt_weights[i]*sqrt_weights[j];
          }
        }

        return CovMat;
      }

    We see that for the diagonal elements of the covariance matrix, we must add the squared
    statistical error. This is done by creating a diagonal matrix from the vector of statistical
    errors, called ``stat_mat``, and squaring the matrix before the return line.

    We now pick datapoint i and loop over all other datapoints j. For each such pair of points,
    loop over all the systematics, if the l'th systematic of data point i has name ``SKIP`` we ignore
    it. This is handled by scanning over all sytematic errors in the ``sys_errors`` dataframe and dropping
    any columns which correspond to a systematic error name of ``SKIP``, thus the ``sys_errors`` dataframe
    defined below contains only the systematics without the ``SKIP`` name.

    Note that in the switch statement an ADDitive or MULTiplicative systype of systematic i is handled
    by either multiplying the additive or multiplicative uncertainties respectively. We instead opt to
    create a purely additive systematic table where all elements are to be understood as additive
    systematics and systematics of type MULT have been correctly transformed to their additive
    values. This dataframe we call ``additive_mat``.

    Next we notice that a systematic is correlated if its name is neither ``UNCORR`` nor ``THEORYUNCORR``. We
    create a dataframe called ``correlated_mat`` where any column corresponding to an ``UNCORR`` or
    ``THEORYUNCORR`` systematic is set to 0.

    From here it's a matter of staring at a piece of paper for a while to realise the contribution to the
    covariance matrix arising due to correlated systematics is ``correlated_mat @ additive_mat.T``. Note
    that in the case where ``i == j`` is met, we double count this contribution using this method, so we set
    the diagonal elements of this particular contribution to 0 by using ``np.fill_diagonal``.

    Finally, the ``i == j`` case is handled by summing the square row elements of ``additive_mat``, this is
    done by using the ``np.einsum`` function.

    For more information on the generation of the covariance matrix see the `paper <https://arxiv.org/pdf/hep-ph/0501067.pdf>`_
    outlining the procedure, specifically equation 2 and surrounding text.

    Paramaters
    ----------
    commondata : validphys.coredata.CommonData
        CommonData which stores information about systematic errors,
        their treatment and description.

    Returns
    -------
    cov_mat : np.array
        Numpy array which is N_dat x N_dat (where N_dat is the number of data points after cuts)
        containing uncertainty and correlation information.

    Example
    -------
    >>> from validphys.commondataparser import load_commondata
    >>> from validphys.loader import Loader
    >>> from validphys.calcutils import covmat_from_systematics
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC")
    >>> cd = load_commondata(cd)
    >>> covmat_from_systematics(cd)
    array([[8.64031971e-05, 8.19971921e-05, 6.27396915e-05, ...,
            2.40747732e-05, 2.79614418e-05, 3.46727332e-05],
           [8.19971921e-05, 1.41907442e-04, 6.52360141e-05, ...,
            2.36624379e-05, 2.72605623e-05, 3.45492831e-05],
           [6.27396915e-05, 6.52360141e-05, 9.41928691e-05, ...,
            1.79244824e-05, 2.08603130e-05, 2.56283708e-05],
           ...,
           [2.40747732e-05, 2.36624379e-05, 1.79244824e-05, ...,
            5.67822050e-05, 4.09077450e-05, 4.14126235e-05],
           [2.79614418e-05, 2.72605623e-05, 2.08603130e-05, ...,
            4.09077450e-05, 5.55150870e-05, 4.15843357e-05],
           [3.46727332e-05, 3.45492831e-05, 2.56283708e-05, ...,
            4.14126235e-05, 4.15843357e-05, 1.43824457e-04]])
    """
    systype = commondata.systype_table["name"]
    # Dropping systematics that have type SKIP
    sys_errors = commondata.sys_errors.loc[:, systype != "SKIP"]

    # Diagonal matrix containing the statistical uncertainty for each
    # data point
    stat_mat = np.diag(commondata.stat_errors.to_numpy())

    # Systematic uncertainties converted to additive uncertainties
    additive_mat = sys_errors.apply(
        lambda x: [
            i.add if i.sys_type == "ADD" else (i.mult * j / 100)
            for i, j in zip(x, commondata.central_values)
        ]
    )

    # Matrix where all non-zero values are due to correlated systematics
    correlated_mat = additive_mat.copy()
    # If all systypes were SKIP then we'd have an empty df
    if not correlated_mat.empty:
        correlated_mat.loc[:, (systype == "UNCORR") | (systype == "THEORYUNCORR")] = 0
        correlated_mat = correlated_mat.to_numpy()
        additive_mat = additive_mat.to_numpy()
    else:
        # Some datasets have zero systematic uncertainties.
        # We treat this special edge case by creating
        # a matrix of zeros so it can be broadcast
        # with the statistical error matrix at the end
        correlated_mat = additive_mat = np.zeros_like(stat_mat)

    # Avoid double counting in case that the i = j element is a correlated
    # systematic. We take care of this case in the einsum below
    off_diagonals = correlated_mat @ additive_mat.T
    np.fill_diagonal(off_diagonals, 0)

    cov_mat = (
        stat_mat ** 2
        + off_diagonals
        + np.diag(np.einsum("ij, ij -> i", additive_mat, additive_mat))
    )
    return cov_mat


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
    >>> from validphys.commondataparser import load_commondata
    >>> from validphys.loader import Loader
    >>> from validphys.calcutils import covmat_from_systematics
    >>> from validphys.results import sqrt_covmat
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC")
    >>> cd = load_commondata(cd)
    >>> cov = covmat_from_systematics(cd)
    >>> sqrtcov = sqrt_covmat(cov)
    array([[9.29533200e-03, 0.00000000e+00, 0.00000000e+00, ...,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
           [8.82133011e-03, 8.00572153e-03, 0.00000000e+00, ...,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
           [6.74959124e-03, 7.11446377e-04, 6.93755946e-03, ...,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
           ...,
           [2.58998530e-03, 1.01842488e-04, 5.34315873e-05, ...,
            4.36182637e-03, 0.00000000e+00, 0.00000000e+00],
           [3.00811652e-03, 9.05569142e-05, 7.09658356e-05, ...,
            1.18572366e-03, 4.31943367e-03, 0.00000000e+00],
           [3.73012316e-03, 2.05432491e-04, 4.40226875e-05, ...,
            9.61421910e-04, 5.31447414e-04, 9.98677667e-03]])
    """
    dimensions = covariance_matrix.shape

    if covariance_matrix.size == 0:
        raise ValueError("Attempting the decomposition of an empty matrix.")
    elif dimensions[0] != dimensions[1]:
        raise ValueError("The input covariance matrix should be square but "
                         f"instead it has dimensions {dimensions[0]} x "
                         f"{dimensions[1]}")

    sqrt_diags = np.sqrt(np.diag(covariance_matrix))
    correlation_matrix = covariance_matrix / sqrt_diags[:, np.newaxis] / sqrt_diags
    decomp = la.cholesky(correlation_matrix)
    sqrt_matrix = (decomp * sqrt_diags).T
    return sqrt_matrix


def groups_covmat_no_table(
       groups_data, groups_index, groups_covmat_collection):
    """Export the covariance matrix for the groups. It exports the full
    (symmetric) matrix, with the 3 first rows and columns being:

        - group name

        - dataset name

        - index of the point within the dataset.
    """
    data = np.zeros((len(groups_index),len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_covmat in zip(
            groups_data, groups_covmat_collection):
        name = group.name
        df.loc[[name],[name]] = group_covmat
    return df


@table
def groups_covmat(groups_covmat_no_table):
    """Duplicate of groups_covmat_no_table but with a table decorator."""
    return groups_covmat_no_table


@table
def groups_sqrtcovmat(
        groups_data, groups_index, groups_sqrt_covmat):
    """Like groups_covmat, but dump the lower triangular part of the
    Cholesky decomposition as used in the fit. The upper part indices are set
    to zero.
    """
    data = np.zeros((len(groups_index),len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_sqrt_covmat in zip(
            groups_data, groups_sqrt_covmat):
        name = group.name
        group_sqrt_covmat[np.triu_indices_from(group_sqrt_covmat, k=1)] = 0
        df.loc[[name],[name]] = group_sqrt_covmat
    return df


@table
def groups_invcovmat(
        groups_data, groups_index, groups_covmat_collection):
    """Compute and export the inverse covariance matrix.
    Note that this inverts the matrices with the LU method which is
    suboptimal."""
    data = np.zeros((len(groups_index),len(groups_index)))
    df = pd.DataFrame(data, index=groups_index, columns=groups_index)
    for group, group_covmat in zip(
            groups_data, groups_covmat_collection):
        name = group.name
        #Improve this inversion if this method tuns out to be important
        invcov = la.inv(group_covmat)
        df.loc[[name],[name]] = invcov
    return df


@table
def groups_normcovmat(groups_covmat, groups_data_values):
    """Calculates the grouped experimental covariance matrix normalised to data."""
    df = groups_covmat
    groups_data_array = np.array(groups_data_values)
    mat = df/np.outer(groups_data_array, groups_data_array)
    return mat


@table
def groups_corrmat(groups_covmat):
    """Generates the grouped experimental correlation matrix with groups_covmat as input"""
    df = groups_covmat
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
    return mat


@check_data_cuts_match_theorycovmat
def dataset_inputs_covmat(
        data: DataGroupSpec,
        fitthcovmat,
        t0set:(PDF, type(None)) = None,
        norm_threshold=None):
    """Like `covmat` except for a group of datasets"""
    loaded_data = data.load()

    if t0set:
        #Copy data to avoid chaos
        loaded_data = type(loaded_data)(loaded_data)
        log.debug("Setting T0 predictions for %s" % data)
        loaded_data.SetT0(t0set.load_t0())

    covmat = loaded_data.get_covmat()

    if fitthcovmat:
        loaded_thcov = fitthcovmat.load()
        ds_names = loaded_thcov.index.get_level_values(1)
        indices = np.in1d(ds_names, [ds.name for ds in data.datasets]).nonzero()[0]
        covmat += loaded_thcov.iloc[indices, indices].values
    if norm_threshold is not None:
        covmat = regularize_covmat(
            covmat,
            norm_threshold=norm_threshold
        )
    return covmat


@check_dataset_cuts_match_theorycovmat
def covmat(
    dataset:DataSetSpec,
    fitthcovmat,
    t0set:(PDF, type(None)) = None,
    norm_threshold=None):
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
    `calcutils.regularize_covmat` if the user specifies `norm_threshold. This
    algorithm sets a minimum threshold for eigenvalues that the corresponding
    correlation matrix can have to be:

    1/(norm_threshold)^2

    which has the effect of limiting the L2 norm of the inverse of the correlation
    matrix. By default norm_threshold is None, to which means no regularization
    is performed.

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
    covmat : array
        a covariance matrix as a numpy array

    Examples
    --------

    >>> from validphys.api import API
    >>> inp = {
            'dataset_input': {'ATLASTTBARTOT'},
            'theoryid': 52,
            'use_cuts': 'no_cuts'
        }
    >>> cov = API.covmat(**inp)
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
    if norm_threshold is not None:
        covmat = regularize_covmat(
            covmat,
            norm_threshold=norm_threshold
        )
    return covmat


@check_pdf_is_montecarlo
def pdferr_plus_covmat(dataset, pdf, covmat):
    """For a given `dataset`, returns the sum of the covariance matrix given by
    `covmat` and the PDF error: a covariance matrix estimated from the
    replica theory predictions from a given monte carlo `pdf`

    Parameters
    ----------
    dataset: DataSetSpec
        object parsed from the `dataset_input` runcard key
    pdf: PDF
        monte carlo pdf used to estimate PDF error
    covmat: np.array
        experimental covariance matrix

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

    See Also
    --------
    covmat: Standard experimental covariance matrix
    """
    loaded_data = dataset.load()
    th = ThPredictionsResult.from_convolution(pdf, dataset, loaded_data=loaded_data)
    pdf_cov = np.cov(th._rawdata, rowvar=True)
    return pdf_cov + covmat


def pdferr_plus_dataset_inputs_covmat(data, pdf, dataset_inputs_covmat):
    """Like `pdferr_plus_covmat` except for an experiment"""
    # do checks get performed here?
    return pdferr_plus_covmat(data, pdf, dataset_inputs_covmat)


def dataset_inputs_sqrt_covmat(dataset_inputs_covariance_matrix):
    """Like `sqrt_covmat` but for an group of datasets"""
    return sqrt_covmat(dataset_inputs_covariance_matrix)


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
    each_dataset, datasets_covmat_no_reg, datasets_covmat_reg, norm_threshold):
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


@check_speclabels_different
@table
def dataspecs_datasets_covmat_differences_table(
    dataspecs_speclabel, dataspecs_covmat_diff_tables
):
    """For each dataspec calculate and tabulate the two covmat differences
    described in `datasets_covmat_differences_table`
    (max relative difference in variance and max absolute correlation difference)

    """
    df = pd.concat( dataspecs_covmat_diff_tables, axis=1)
    cols = df.columns.get_level_values(0).unique()
    df.columns = pd.MultiIndex.from_product((dataspecs_speclabel, cols))
    return df


groups_covmat_collection = collect(
    'dataset_inputs_covariance_matrix', ('group_dataset_inputs_by_metadata',)
)

groups_sqrt_covmat = collect(
    'dataset_inputs_sqrt_covmat',
    ('group_dataset_inputs_by_metadata',)
)

dataspecs_covmat_diff_tables = collect(
    "datasets_covmat_differences_table", ("dataspecs",)
)

fits_name_with_covmat_label = collect('fit_name_with_covmat_label', ('fits',))

datasets_covmat_no_reg = collect(
    "covariance_matrix", ("data", "no_covmat_reg"))

datasets_covmat_reg = collect(
    "covariance_matrix", ("data",))

datasets_covmat = collect('covariance_matrix', ('data',))

datasets_covariance_matrix = collect(
    'covariance_matrix',
    ('experiments', 'experiment',)
)
