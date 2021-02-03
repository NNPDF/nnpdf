"""
covmat_construction.py

Module containing underlying functions which construct the covariance matrices

"""
import numpy as np
import pandas as pd
import scipy.linalg as la

INTRA_DATASET_SYS_NAME = ("UNCORR", "CORR", "THEORYUNCORR", "THEORYCORR")


def construct_covmat(stat_errors: np.array, sys_errors: pd.DataFrame):
    """Basic function to construct a covariance matrix (covmat), given the
    statistical error and a dataframe of systematics.

    Errors with name UNCORR or THEORYUNCORR are added in quadrature with
    the statistical error to the diagonal of the covmat.

    Other systematics are treated as correlated; their covmat contribution is
    found by multiplying them by their transpose.

    Parameters
    ----------
    stat_errors: np.array
        a 1-D array of statistical uncertainties
    sys_errors: pd.DataFrame
        a dataframe with shape (N_data * N_sys) and systematic name as the
        column headers. The uncertainties should be in the same units as the
        data.

    Notes
    -----
    This function doesn't contain any logic to ignore certain contributions to
    the covmat, if you wanted to not include a particular systematic/set of
    systematics i.e all uncertainties with MULT errors, then filter those out
    of ``sys_errors`` before passing that to this function.

    """
    diagonal = stat_errors ** 2

    is_uncorr = sys_errors.columns.isin(("UNCORR", "THEORYUNCORR"))
    diagonal += (sys_errors.loc[:, is_uncorr].to_numpy() ** 2).sum(axis=1)

    corr_sys_mat = sys_errors.loc[:, ~is_uncorr].to_numpy()
    return np.diag(diagonal) + corr_sys_mat @ corr_sys_mat.T


def covmat_from_systematics(commondata, central_values=None):
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
    commondata : validphys.coredata.CommonData
        CommonData which stores information about systematic errors,
        their treatment and description.
    central_values : None, np.array
        1-D array containing alternative central values to combine with the
        multiplicative errors to calculate their absolute contributions. By
        default this is None, and the experimental central values are used. However, this
        can be used to calculate, for example, the t0 covariance matrix by
        using the predictions from the central member of the t0 pdf.

    Returns
    -------
    cov_mat : np.array
        Numpy array which is N_dat x N_dat (where N_dat is the number of data
        points after cuts) containing uncertainty and correlation information.

    Example
    -------
    >>> from validphys.results_providers.commondata_parser.py import load_commondata
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
    return construct_covmat(
        commondata.stat_errors.to_numpy(),
        commondata.systematic_errors(central_values)
    )


def datasets_covmat_from_systematics(
    list_of_commondata, list_of_central_values=None
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
    list_of_commondata : list[validphys.coredata.CommonData]
        list of CommonData objects.
    list_of_central_values: None, list[np.array]
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
    >>> from validphys.results_providers.commondata_parser.py import load_commondata
    >>> from validphys.covmats import datasets_covmat_from_systematics
    >>> from validphys.loader import Loader
    >>> l = Loader()
    >>> cd1 = l.check_commondata("ATLASLOMASSDY11EXT")
    >>> cd2 = l.check_commondata("ATLASZHIGHMASS49FB")
    >>> ld1, ld2 = map(load_commondata, (cd1, cd2))
    >>> datasets_covmat_from_systematics((ld1, ld2))
    array([[2.91814548e+06, 4.66692123e+06, 2.36823008e+06, 8.62587330e+05,
            2.78209614e+05, 1.11790645e+05, 1.75129920e+03, 7.97466600e+02,
            4.00296960e+02, 2.22039720e+02, 1.46202210e+02, 8.36558100e+01,
    """
    special_corrs = []
    block_diags = []

    if list_of_central_values is None:
        # want to just pass None to systematic_errors method
        list_of_central_values = [None] * len(list_of_commondata)

    for cd, central_values in zip(list_of_commondata, list_of_central_values):
        errors = cd.systematic_errors(central_values)
        # separate out the special uncertainties which can be correlated across
        # datasets
        is_intra_dataset_error = errors.columns.isin(INTRA_DATASET_SYS_NAME)
        block_diags.append(construct_covmat(
            cd.stat_errors.to_numpy(), errors.loc[:, is_intra_dataset_error]))
        special_corrs.append(errors.loc[:, ~is_intra_dataset_error])

    # concat systematics across datasets
    special_sys = pd.concat(special_corrs, axis=0, sort=False)
    # non-overlapping systematics are set to NaN by concat, fill with 0 instead.
    special_sys.fillna(0, inplace=True)

    diag = la.block_diag(*block_diags)
    return diag + special_sys.to_numpy() @ special_sys.to_numpy().T
