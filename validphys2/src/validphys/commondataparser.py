"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
from operator import attrgetter

import numpy as np
import pandas as pd
import scipy.linalg as la

from validphys.core import peek_commondata_metadata
from validphys.coredata import CommonData, SystematicError

def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    Returns an instance of CommonData
    """
    commondatafile = spec.datafile
    # Getting set name from commondata file name
    setname = commondatafile.name[5:-4] # DATA prefix and .dat suffix
    systypefile = spec.sysfile

    commondata = parse_commondata(commondatafile, systypefile, setname)

    return commondata


def parse_commondata(commondatafile, systypefile, setname):
    """Parse a commondata file  and a systype file into a CommonData.

    Parameters
    ----------
    commondatafile : file or path to file
    systypefile : file or path to file

    Returns
    -------
    commondata : CommonData
        An object containing the data and information from the commondata
        and systype files.
    """
    # First parse commondata file
    commondatatable = pd.read_csv(commondatafile, sep=r'\s+', skiprows=1, header=None)
    # Remove NaNs
    # TODO: replace commondata files with bad formatting
    # Build header
    commondataheader = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    nsys  = (commondatatable.shape[1] - len(commondataheader)) // 2
    for i in range(nsys):
        commondataheader += [f"sys.add.{i+1}", f"sys.mult.{i+1}"]
    commondatatable.columns = commondataheader
    commondatatable.set_index("entry", inplace=True)
    ndata = len(commondatatable)
    commondataproc = commondatatable["process"][1]
     # Check for consistency with commondata metadata
    cdmetadata =  peek_commondata_metadata(commondatafile)
    if (setname, nsys, ndata) != attrgetter('name', 'nsys', 'ndata')(cdmetadata):
        raise ValueError("Commondata table information does not match metadata")

    # Now parse the systype file
    systypetable = parse_systypes(systypefile)

    # Populate CommonData object
    return CommonData(
        setname=setname,
        ndata=ndata,
        commondataproc=commondataproc,
        nkin=3,
        nsys=nsys,
        commondata_table=commondatatable,
        systype_table=systypetable
    )

def parse_systypes(systypefile):
    """Parses a systype file and returns a pandas dataframe.
    """
    systypeheader = ["sys_index", "type", "name"]
    try:
        systypetable = pd.read_csv(
            systypefile, sep=r"\s+", names=systypeheader, skiprows=1, header=None
        )
        systypetable.dropna(axis='columns', inplace=True)
    # Some datasets e.g. CMSWCHARMRAT have no systematics
    except pd.errors.EmptyDataError:
        systypetable = pd.DataFrame(columns=systypeheader)

    systypetable.set_index("sys_index", inplace=True)

    return systypetable


def combine_commondata(commondata_list):
    """Function that takes in a sequence of :py:class:`validphys.coredata.CommonData`
    objects and returns an effective :py:class:`validphys.coredata.CommonData` object
    which has taken into account the (possible) correlations between the datasets. This
    effective ``CommonData`` can then be used to obtain a correlated covariance matrix,
    for example.

    .. note::
        The return of this function does not have a file on disk associated with it.
        As such the ``CommonData`` instance has a ``setname`` of ``"N/A"`` and a
        ``processtype`` of ``"MIXED"``.

    Parameters
    ----------
    commondata_list: sequence
        Sequence of :py:class:`validphys.coredata.CommonData` objects

    Returns
    -------
    comb_cd: validphys.coredata.CommonData
        The combined ``CommonData``

    Raises
    ------
    ValueError: if one or more of the input ``CommonData`` objects have a different
        number of kinematic variables

    Example
    -------
    >>> from validphys.commondataparser import combine_commondata, load_commondata
    >>> from validphys.loader import Loader
    >>> l = Loader()
    >>> cd1 = l.check_commondata("ATLASLOMASSDY11EXT")
    >>> cd2 = l.check_commondata("ATLASZHIGHMASS49FB")
    >>> ld1, ld2 = map(load_commondata, (cd1, cd2))
    >>> comb_cd = combine_commondata((ld1, ld2))
    # Note that the two datasets have one pair of correlated systematics
    # so the total number of systematics is one fewer than the sum
    >>> comb_cd.nsys
    18
    >>> ld1.nsys + ld2.nsys
    19
    >>> from validphys.covmats import covmat_from_systematics
    >>> covmat = covmat_from_systematics(comb_cd)
    >>> covmat.shape
    (19, 19)
    >>> ld1.ndata + ld2.ndata
    19
    """
    # Raise an exception if any of the CommonDatas has
    # a different number of kinematic variables to the others
    # TODO: handle this case instead of raising an error
    nkins = [cd.nkin for cd in commondata_list]
    if not all(x == nkins[0] for x in nkins):
        raise ValueError("One or more of the datasets have a different number of kinematics.")

    # Combine the systematics from each of the CommonDatas
    systematics, systype_table = combine_commondata_systematics(commondata_list)

    # The number of data points is the sum of data points per CommonData while
    # the systematics have correlations accounted for
    ndata, nsys = systematics.shape

    header = ['process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    kinematics = pd.concat([i.commondata_table.loc[:, header] for i in commondata_list], ignore_index=True)
    kinematics.index += 1

    additive = systematics.applymap(lambda x: x.add)
    multiplicative = systematics.applymap(lambda x: x.mult)
    additive.columns = (f"sys.add.{i + 1}" for i in range(nsys))
    multiplicative.columns = (f"sys.mult.{i + 1}" for i in range(nsys))

    # Table containing alternating additive and multiplicative systematic uncertainties
    systematics_header = np.empty((additive.columns.size + multiplicative.columns.size,), dtype=object)
    systematics_header[0::2], systematics_header[1::2] = additive.columns, multiplicative.columns
    systematics_table = additive.join(multiplicative)[systematics_header]

    commondata_table = pd.concat([kinematics, systematics_table], axis=1)

    # The effective CommonData object to return. The set has no unique name or
    # process type and so there is no file associated with this object.
    comb_cd = CommonData("N/A", ndata, "MIXED", 3, nsys, commondata_table, systype_table)

    return comb_cd


def combine_commondata_systematics(commondata_list):
    """Function that takes in a sequence of :py:class:`validphys.coredata.CommonData`
    and returns a tuple of ``pd.DataFrame``. The first ``DataFrame`` contains
    :py:class:`validphys.coredata.SystematicError` which has accounted for the correlations
    between systematics while the second contains information about the systematic errors
    themselves such as their type and name.

    Parameters
    ----------
    commondata_list: sequence
        Sequence of :py:class:`validphys.coredata.CommonData` objects

    Returns
    -------
    eff_sys_errors, systype_table: (pd.DataFrame, pd.DataFrame)
        Tuple of DataFrames containing systematic errors and information about their name and type
    """
    special_types = frozenset({"UNCORR", "CORR", "THEORYUNCORR", "THEORYCORR", "SKIP"})

    # Create an effective systematic error table which is block diagonal
    # where the off-diagonal blocks are filled with zeros as sentinel values.
    eff_sys_errors = pd.DataFrame(la.block_diag(*[i.sys_errors for i in commondata_list]))
    # Replace 0s with np.nan so we can use df.fillna later on
    eff_sys_errors.replace(0, np.nan, inplace=True)
    # Make index and columns start from 1
    eff_sys_errors.index += 1
    eff_sys_errors.columns += 1

    systypes = pd.concat([cd.systype_table for cd in commondata_list], ignore_index=True)
    systypes.index += 1

    # If the systypes are duplicated and not one of the special
    # types above, then those set of systematics are correlated
    # and we keep the first occurence and replace each subsequent
    # occurence with a sentinel value of None
    duplicates = systypes.groupby("name").groups
    # Remove the duplicates that are special types
    [duplicates.pop(key, None) for key in special_types]

    for cols in duplicates.values():
        # The columns that were correlated
        to_combine = eff_sys_errors.loc[:, cols]
        # We backwards propagate the uncertainties
        to_combine.fillna(method="bfill", axis=1, inplace=True)
        # Replace the eff_sys_errors column with
        # the systematics that were correlated
        prepared_column = to_combine.iloc[:, 0]
        eff_sys_errors.loc[:, cols[0]] = prepared_column

    # The zeros are replaced with validphys.coredata.SystematicError types, we label the
    # name as "SKIP" so any future logic knows that these systematic errors should be ignored
    dummy_error = SystematicError(0, 0, None, "SKIP")
    eff_sys_errors.fillna(dummy_error, inplace=True)

    keep_columns = (systypes["name"].isin(special_types)) | ~(
        systypes.duplicated(subset="name", keep="first")
    )
    keep_columns.index = eff_sys_errors.columns
    eff_sys_errors = eff_sys_errors.loc[:, keep_columns]
    # After slicing with keep_columns, the index
    # of the systematic errors need to be reset
    eff_sys_errors.columns = range(1, len(eff_sys_errors.columns) + 1)

    systype_table = systypes[keep_columns]
    systype_table.index = range(1, len(eff_sys_errors.columns) + 1)

    return eff_sys_errors, systype_table
