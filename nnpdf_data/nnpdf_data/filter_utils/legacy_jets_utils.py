"""
Utils specifically needed to ease the implementation of the 7 and 8 TeV jet datasets,
these are: CMS_2JET_7TEV, CMS_1JET_8TEV, ATLAS_1JET_8TEV_R06, ATLAS_2JET_7TEV_R06
"""

import numpy as np
import pandas as pd
from scipy.linalg import block_diag
import yaml

def range_str_to_floats(str_range):
    """
    converts a string range to a list,
    e.g. "0.5 - 1.0" --> [0.5,1.0]
    and returns a dictionary
    """
    # Split the range string into two parts
    str_nums = str_range.split('-')
    # Convert the parts to floats
    min = float(str_nums[0])
    max = float(str_nums[1])
    mid = float(f"{0.5 * (min + max):.3f}")
    # Return a dict
    return {"min": min, "mid": mid, "max": max}


# ==================================================================== CMS_2JET_7TEV ====================================================================#

def get_data_values_CMS_2JET_7TEV(tables, version):
    """
    returns the central data

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []
    for table in tables:

        hepdata_table = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']

        for value in values:
            data_central.append(value['value'])

    return data_central


def get_kinematics_CMS_2JET_7TEV(tables, version):
    """
    returns the relevant kinematics values

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    for table in tables:

        hepdata_table = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        rapidity_interval = input['dependent_variables'][0]['qualifiers'][0]['value']
        ydiff = {}
        if rapidity_interval == '< 0.5':
            ydiff['min'], ydiff['max'], ydiff['mid'] = 0.0, 0.5, 0.25
        else:
            ydiff = range_str_to_floats(rapidity_interval)

        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])

        dijet_mass_bins = input['independent_variables'][0]['values']

        for m_jj in dijet_mass_bins:

            # kinematics
            m_jj['low'], m_jj['high'] = m_jj['low'], m_jj['high']
            m_jj['mid'] = float(f"{0.5 * (m_jj['low']+m_jj['high']):.3f}")

            kin_value = {
                'ydiff': {'min': ydiff['min'], 'mid': ydiff['mid'], 'max': ydiff['max']},
                'm_jj': {'min': m_jj['low'], 'mid': m_jj['mid'], 'max': m_jj['high']},
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)

    return kin


def get_corr_dat_file_CMS_2JET_7TEV(filename):
    """
    read out correlation matrices from the
    dijet_corr.dat file

    Parameters
    ----------
    filename : str
        Takes path to dijet_corr.dat

    Returns
    -------
    list
        list, each element of which is a 2D np array
    """

    with open(filename) as file:
        lines = file.readlines()

    # store the number of the rows where the correlation matrix
    # is printed
    begin_rows = []
    end_rows = []

    for i, line in enumerate(lines):

        if "Statistical correlation" in line and begin_rows == []:
            begin_rows.append(i + 2)

        elif "Statistical correlation" in line:
            begin_rows.append(i + 2)
            end_rows.append(i - 2)

        elif i == len(lines) - 1:
            end_rows.append(i)

    correlation_matrices = []
    for begin_row, end_row in zip(begin_rows, end_rows):

        size_mat = end_row - begin_row + 1
        stat_corr = np.zeros((size_mat, size_mat))

        i = 0
        for idx in range(begin_row, end_row + 1):
            # ignore first two columns as these give the bin kin
            stat_corr[i] = np.fromstring(lines[idx], sep=' ')[2:]
            i += 1

        correlation_matrices.append(stat_corr)

    return correlation_matrices


def get_stat_uncertainties_CMS_2JET_7TEV():
    """
    function used to get the statistical
    uncertainty from the HEPdata tables.

    Returns
    -------
    dict
        dictionary with keys = number of table
        value = list of statistical uncertainties

    """

    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    stat_err = {}

    for table in tables:
        stat = []
        hepdata_tables = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"
        with open(hepdata_tables) as file:
            input = yaml.safe_load(file)

        for err in input['dependent_variables'][0]['values']:
            stat.append(err['errors'][0]['symerror'])
        stat_err[table] = stat

    return stat_err


def dat_file_to_df_CMS_2JET_7TEV():
    """
    from dijet_sys.dat table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and
    np.nan entries

    Returns
    -------
    list
        list of dataframes

    """

    with open("rawdata/dijet_sys.dat") as file:
        lines = file.readlines()

    # get rows of numeric tables in dat file
    begin_lines = []
    end_lines = []

    for i, line in enumerate(lines):

        if "mlo" in line:
            begin_lines.append(i + 2)

        if "The NP correction" in line and len(begin_lines) >= 1:
            end_lines.append(i - 1)

    end_lines.append(len(lines))

    # define dataframe
    columns = [
        "JEC0-",
        "JEC0+",
        "JEC1-",
        "JEC1+",
        "JEC2-",
        "JEC2+",
        "JEC3-",
        "JEC3+",
        "JEC4-",
        "JEC4+",
        "JEC5-",
        "JEC5+",
        "JEC6-",
        "JEC6+",
        "JEC7-",
        "JEC7+",
        "JEC8-",
        "JEC8+",
        "JEC9-",
        "JEC9+",
        "JEC10-",
        "JEC10+",
        "JEC11-",
        "JEC11+",
        "JEC12-",
        "JEC12+",
        "JEC13-",
        "JEC13+",
        "Lumi-",
        "Lumi+",
        "Unfolding+",
        "Unfolding-",
        "Bin-by-bin+",
        "Bin-by-bin-",
    ]

    dfs = []
    for beg, end in zip(begin_lines, end_lines):

        df = pd.DataFrame(columns=columns, index=range(beg, end))
        j = 0
        for i in range(beg, end):
            # do not consider NP uncertainty
            col_vals = np.fromstring(lines[i], sep=' ')[5:]

            df.iloc[j] = col_vals
            j += 1

        dfs.append(df)

    return dfs


def JEC_error_matrix_CMS_2JET_7TEV():
    """
    Jet Energy Scale (JET): 14 Asymmetric uncertainties correlated across all
    rapidity and mass bins (CORR). This uncertainty is not always presented as
    (left<0 and right>0), e.g. [-delta_left,+delta_right]
    Hence the D'Agostini prescription for symmetrising errors
    is not valid here because it works with the only case displayed above.
    Instead, we use here the experimentalists prescription, where we take every
    subpart of the uncertainty left and right as independent source of
    uncertainty. This is motivated by taking the average
    of the left and right uncertainty, hence the origin of the sqrt(2)
    that we divide by.


    Returns
    -------
    np.array
        covariance matrix for JET energy scale uncertainty

    """
    dfs = dat_file_to_df_CMS_2JET_7TEV()
    JEC_err = []
    for df in dfs:
        JEC_err.append(df.filter(like="JEC"))

    # divide by sqrt(2) since treating each unc of asymm as independent
    jec = pd.concat(JEC_err, axis=0) / np.sqrt(2)

    # get central value to convert mult error
    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']
    cv = get_data_values_CMS_2JET_7TEV(tables, version)
    cv = np.array(cv)

    # convert mult error to absolute
    jec = jec.multiply(cv, axis=0)

    return jec


def lumi_covmat_CMS_2JET_7TEV():
    """
    Luminosity uncertainty: this is a symmetric uncertainty of 2.2% correlated
    accross all mass and rapidity bins and all CMS datasets at 7 TeV, hence the
    keyword (CMSLUMI11).

    NOTE: this function is needed to test only whether the full covmat coincides
    with the old implementation and can be removed at some point

    Returns
    -------
    np.array
        covariance matrix for luminosity uncertainty

    """
    dfs = dat_file_to_df_CMS_2JET_7TEV()
    lumi_err = []
    for df in dfs:
        lumi_err.append(df.filter(like="Lumi+"))

    lumi = pd.concat(lumi_err, axis=0)

    # get central value to convert mult error
    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']
    cv = get_data_values_CMS_2JET_7TEV(tables, version)
    cv = np.array(cv)

    # convert mult to abs
    lumi = lumi.multiply(cv, axis=0)
    lumi = lumi.to_numpy().astype(float)

    return np.einsum('ij,kj->ik', lumi, lumi)


def unfolding_error_matrix_CMS_2JET_7TEV():
    """
    Unfolding uncertainty: this asymmetric is correlated across all rapidity
    and mass bins (CORR).

    Returns
    -------
    np.array
        covariance matrix for unfolding uncertainty
    """
    dfs = dat_file_to_df_CMS_2JET_7TEV()
    unfold_err = []
    for df in dfs:
        unfold_err.append(df.filter(like="Unfolding"))

    unfold = pd.concat(unfold_err, axis=0) / np.sqrt(2)

    # get central value to convert mult error
    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']
    cv = get_data_values_CMS_2JET_7TEV(tables, version)
    cv = np.array(cv)

    # convert mult to abs
    unfold = unfold.multiply(cv, axis=0)

    return unfold


def bin_by_bin_covmat_CMS_2JET_7TEV():
    """
    Bin-by-Bin uncertainty: this is a symmetric uncertainty fully uncorrelated
    accross bins of mass and rapidity (UNCORR)

    NOTE: this function is needed to test only whether the full covmat coincides
    with the old implementation and can be removed at some point

    Returns
    -------
    np.array
        covariance matrix for bin by bin uncertainty
    """
    dfs = dat_file_to_df_CMS_2JET_7TEV()
    bin_err = []
    for df in dfs:
        bin_err.append(df.filter(like="Bin-by-bin-"))  # symm so choose only one

    bin = pd.concat(bin_err, axis=0)

    # get central value to convert mult error
    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']
    cv = get_data_values_CMS_2JET_7TEV(tables, version)
    cv = np.array(cv)

    # convert mult to abs
    bin = bin.multiply(cv, axis=0)

    bin = bin.to_numpy().astype(float)

    # fully uncorrelated
    bin_cov = np.diag(bin.reshape(bin.shape[0]) ** 2)
    return bin_cov


# ==================================================================== CMS_1JET_8TEV ====================================================================#


TABLE_TO_RAPIDITY_CMS_1JET_8TEV = {
    1: [0.0, 0.5],
    2: [0.5, 1.0],
    3: [1.0, 1.5],
    4: [1.5, 2],
    5: [2.0, 2.5],
    6: [2.5, 3.0],
}

COLUMN_NAMES_CMS_1JET_8TEV = [
    'binFlag',
    'ylow',
    'yhigh',
    'ptlow',
    'pthigh',
    'Sigma',
    'NPCorr',
    'npcorerr+',
    'npcorerr-',
    'ignore',
    'Unfolding+',
    'Unfolding-',
    'AbsoluteStat+',
    'AbsoluteStat-',
    'AbsoluteScale+',
    'AbsoluteScale-',
    'AbsoluteMPFBias+',
    'AbsoluteMPFBias-',
    'Fragmentation+',
    'Fragmentation-',
    'SinglePionECAL+',
    'SinglePionECAL-',
    'SinglePionHCAL+',
    'SinglePionHCAL-',
    'FlavorQCD+',
    'FlavorQCD-',
    'RelativeJEREC1+',
    'RelativeJEREC1-',
    'RelativeJEREC2+',
    'RelativeJEREC2-',
    'RelativeJERHF+',
    'RelativeJERHF-',
    'RelativePtBB+',
    'RelativePtBB-',
    'RelativePtEC1+',
    'RelativePtEC1-',
    'RelativePtEC2+',
    'RelativePtEC2-',
    'RelativePtHF+',
    'RelativePtHF-',
    'RelativeFSR+',
    'RelativeFSR-',
    'RelativeStatEC2+',
    'RelativeStatEC2-',
    'RelativeStatHF+',
    'RelativeStatHF-',
    'PileUpDataMC+',
    'PileUpDataMC-',
    'PileUpPtRef+',
    'PileUpPtRef-',
    'PileUpPtBB+',
    'PileUpPtBB-',
    'PileUpPtEC1+',
    'PileUpPtEC1-',
    'PileUpPtEC2+',
    'PileUpPtEC2-',
    'PileUpPtHF+',
    'PileUpPtHF-',
    'RelativeStatFSR+',
    'RelativeStatFSR-',
    'Luminosity',
    'stat',
    'uncor',
]

TABLE_DATA_SHAPE_CMS_1JET_8TEV = {1: 37, 2: 37, 3: 36, 4: 32, 5: 25, 6: 18}


def get_data_values_CMS_1JET_8TEV(tables, version):
    """
    returns the central data

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []
    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']

        for value in values:
            data_central.append(value['value'])

    return data_central


def get_kinematics_CMS_1JET_8TEV(tables, version):
    """
    returns the relevant kinematics values

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        # rapidity
        rapidity_interval = TABLE_TO_RAPIDITY_CMS_1JET_8TEV[table]
        rap = {}
        rap['min'], rap['max'] = rapidity_interval[0], rapidity_interval[1]
        rap['mid'] = 0.5 * (rap['min'] + rap['max'])

        # center of mass energy
        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])

        # transverse momentum
        jet_kt_bins = input['independent_variables'][0]['values']
        KT = {}
        for kt in jet_kt_bins:

            KT['min'], KT['max'] = kt['low'], kt['high']
            KT['mid'] = float(f"{0.5 * (kt['low'] + kt['high']):.3f}")

            kin_value = {
                'y': {'min': rap['min'], 'mid': rap['mid'], 'max': rap['max']},
                'pT': {'min': KT['min'], 'mid': KT['mid'], 'max': KT['max']},
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)

    return kin


def get_stat_correlations_CMS_1JET_8TEV(table):
    """
    Read a CMS_8TeV_jets_Ybin#___CMS_8TeV_jets_Ybin#.dat
    file and output the statistical correlation matrix.


    Parameters
    ----------
    table : int
            number of the table

    Returns
    -------
    np.array
            2-D array

    """
    with open(f'rawdata/CMS_8TeV_jets_Ybin{table}___CMS_8TeV_jets_Ybin{table}.dat') as file:
        card = file.readlines()

    # get shape of matrix
    shape_mat = TABLE_DATA_SHAPE_CMS_1JET_8TEV[table]

    stat_corr = np.zeros((shape_mat, shape_mat))

    # correlation rows always start at row 18
    for j in range(shape_mat):
        # fill rows of correlation matrix
        stat_corr[j, :] = np.array(
            [card[(17 + shape_mat * j) + k].split()[-1] for k in range(shape_mat)]
        )

    # add zeros for points in the pt<74 kinematic region
    # these points should be cut (it is always 9 pt bins in the pt < 74 region)
    stat_corr = block_diag(np.zeros((9, 9)), stat_corr)
    return stat_corr


def block_diagonal_corr_CMS_1JET_8TEV(tables):
    """
    forms block diagonal correlation matrix
    for stat uncertainties. Each block corresponds
    to a rapidity bin.

    Parameters
    ----------
    tables : list
            list of integers numbering the tables

    Returns
    -------
    np.array
        block diagonal matrix of dim ndata x ndata
    """
    bd_corr = get_stat_correlations_CMS_1JET_8TEV(tables[0])

    for table in tables[1:]:

        bd_corr = block_diag(bd_corr, get_stat_correlations_CMS_1JET_8TEV(table))

    return bd_corr


def get_stat_uncertainties_CMS_1JET_8TEV():
    """
    function used to get the statistical
    uncertainty from the HEPdata tables.

    Returns
    -------
    np.array
        array with ordered stat errors for all tables

    """

    with open('metadata.yaml') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    stat_err = []

    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"
        with open(hepdata_tables) as file:
            input = yaml.safe_load(file)

        # discard pT < 74 GeV entries

        for err, pt in zip(
            input['dependent_variables'][0]['values'], input['independent_variables'][0]['values']
        ):
            stat_err.append(err['errors'][0]['symerror'])

    return np.array(stat_err)


def get_uncertainties_df_CMS_1JET_8TEV(table):
    """ """

    # read dat file into dataframe by skipping the first 41 metadata rows
    df = pd.read_csv(
        f"rawdata/CMS_8TeV_jets_Ybin{table}.dat",
        sep=r"\s+",
        skiprows=41,
        names=COLUMN_NAMES_CMS_1JET_8TEV,
    )

    # reindex
    df = df.reset_index(drop=True)

    df = df[1:-1]  # discard last one as it is repeated

    return df


def uncertainties_df_CMS_1JET_8TEV(tables):
    """ """
    dfs = []

    for table in tables:
        dfs.append(get_uncertainties_df_CMS_1JET_8TEV(table))
    df = pd.concat(dfs, axis=0)
    return df


def process_err_CMS_1JET_8TEV(df):
    """
    Given the uncertainties dataframe, if the two variations in the pair
    (of uncertainties) have the same sign, only the largest (in absolute value)
    is retained, while the other is set to zero

    """
    for col_idx in np.arange(0, len(df.columns), 2):

        for row_idx, (val1, val2) in enumerate(zip(df.iloc[:, col_idx], df.iloc[:, col_idx + 1])):
            if np.sign(val1) == np.sign(val2):

                if np.abs(val1) > np.abs(val2):
                    df.iloc[row_idx, col_idx + 1] = 0
                elif np.abs(val1) < np.abs(val2):
                    df.iloc[row_idx, col_idx] = 0
    return df


# ============================================================ ATLAS_1JET_8TEV_R06 ============================================================#


TABLE_TO_RAPIDITY_ATLAS_1JET_8TEV_R06 = {
    1: [0.0, 0.5],
    2: [0.5, 1.0],
    3: [1.0, 1.5],
    4: [1.5, 2],
    5: [2.0, 2.5],
    6: [2.5, 3.0],
}


def get_data_values_ATLAS_1JET_8TEV_R06(tables, version):
    """
    returns the central data.
    Note: central data is the same for both correlation scenarios.

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []
    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']

        for value in values:
            data_central.append(value['value'])

    return data_central


def get_kinematics_ATLAS_1JET_8TEV_R06(tables, version):
    """
    returns the relevant kinematics values.
    Note: kinematics are the same for both correlation scenarios.

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        # rapidity
        rapidity_interval = TABLE_TO_RAPIDITY_ATLAS_1JET_8TEV_R06[table]
        rap = {}
        rap['min'], rap['max'] = rapidity_interval[0], rapidity_interval[1]
        rap['mid'] = 0.5 * (rap['min'] + rap['max'])

        # center of mass energy
        sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])

        # transverse momentum
        jet_kt_bins = input['independent_variables'][0]['values']
        KT = {}
        for kt in jet_kt_bins:
            KT['min'], KT['max'] = kt['low'], kt['high']
            KT['mid'] = float(f"{0.5 * (kt['low'] + kt['high']):.3f}")

            kin_value = {
                'y': {'min': rap['min'], 'mid': rap['mid'], 'max': rap['max']},
                'pT': {'min': KT['min'], 'mid': KT['mid'], 'max': KT['max']},
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)

    return kin


def HEP_table_to_df_ATLAS_1JET_8TEV_R06(table, version, variant='nominal'):
    """
    Given hep data table return a pandas DataFrame with index given by Ndata,
    columns by the uncertainties and np.nan entries

    Parameters
    ----------
    table : int
            number of table
    version : int
            version number

    """
    if variant == 'nominal':
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"
    elif variant == 'decorrelated':
        hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"

    with open(hepdata_table) as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][0]['values']

    columns = {}
    errors = card[0]['errors']
    for err in errors:
        # luminosity and stat uncertainty, always symmetric
        if err['label'] == 'syst_lumi' or err['label'] == 'stat':
            columns[err['label']] = np.nan
        else:
            columns[f"{err['label']}_plus"] = np.nan
            columns[f"{err['label']}_minus"] = np.nan

    df_nan = pd.DataFrame(columns, index=range(1, len(card) + 1))

    return df_nan


def fill_df_ATLAS_1JET_8TEV_R06(table, version, variant='nominal'):
    """
    Fill a data frame with index
    corresponding to measured datapoints
    and columns to different uncertainties
    Each df is for a fixed rapidity bin.

    Parameters
    ----------
    table : int
            number of table
    version : int
            version number
    """
    if variant == 'nominal':
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"
    elif variant == 'decorrelated':
        hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"

    with open(hepdata_table) as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][0]['values']
    df = HEP_table_to_df_ATLAS_1JET_8TEV_R06(table, version, variant)

    for i, dat in enumerate(card):
        # cv = dat['value']

        for j, err in enumerate(dat['errors']):
            # lum and stat uncertainties are always symmetric
            if err['label'] == 'stat':
                df.loc[df.index == i + 1, err['label']] = err['asymerror']['plus']
            elif err['label'] == 'syst_lumi':
                df.loc[df.index == i + 1, err['label']] = err['asymerror']['plus']
            else:
                dp, dm = process_err_ATLAS_1JET_8TEV_R06(err)
                df.loc[df.index == i + 1, f"{err['label']}_plus"] = dp  # err['asymerror']['plus']
                df.loc[df.index == i + 1, f"{err['label']}_minus"] = dm  # err['asymerror']['minus']
    return df


def process_err_ATLAS_1JET_8TEV_R06(error):
    """
    Note: the d'Agostini prescription for the
    symmetrization of the error does not hold here.
    We follow the experimentalists prescription.

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}

    Returns
    -------
    tuple
        tuple containing two floats

    """
    d_p = float(error['asymerror']['plus'])
    d_m = float(error['asymerror']['minus'])

    tmp1 = d_p
    tmp2 = d_m
    # case 1: d_p and d_m are both negative
    if tmp1 < 0.0 and tmp2 < 0.0:
        if tmp2 < tmp1:
            d_p = 0.0
            d_m = tmp2
        else:
            d_p = 0.0
            d_m = tmp1

    # case 2: d_p and d_m are both positive
    if tmp1 > 0.0 and tmp2 > 0.0:
        if tmp1 > tmp2:
            d_p = tmp1
            d_m = 0.0
        else:
            d_p = tmp2
            d_m = 0.0
    return d_p, d_m


# ======================================================== ATLAS_2JET_7TEV_R06 ======================================================== #

SCENARIO_ATLAS_2JET_7TEV_R06 = {'nominal': 0, 'stronger': 1, 'weaker': 2}


def process_err_ATLAS_2JET_7TEV_R06(error, cv):
    """
    Converts an error given in percentage
    of the central data value in absolute value.

    Note: the d'Agostini prescription for the
    symmetrization of the error does not hold here.
    We follow here the experimental prescription

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}

    cv : float
        central value

    Returns
    -------
    tuple
        tuple containing two floats

    """
    if "lum" in error:
        # luminosity uncertainty is always symmetric
        sigma = float(error['symerror'].strip('%')) / 100.0 * cv
        return sigma

    elif error['label'] == 'sys':
        if 'asymerror' in error:
            d_p = float(error['asymerror']['plus'].strip('%')) / 100.0 * cv
            d_m = float(error['asymerror']['minus'].strip('%')) / 100.0 * cv

            tmp1 = d_p
            tmp2 = d_m
            # case 1: d_p and d_m are both negative
            if tmp1 < 0.0 and tmp2 < 0.0:
                if tmp2 < tmp1:
                    d_p = 0.0
                    d_m = tmp2
                else:
                    d_p = 0.0
                    d_m = tmp1

            # case 2: d_p and d_m are both positive
            if tmp1 > 0.0 and tmp2 > 0.0:
                if tmp1 > tmp2:
                    d_p = tmp1
                    d_m = 0.0
                else:
                    d_p = tmp2
                    d_m = 0.0
            return d_p / np.sqrt(2.0), d_m / np.sqrt(2.0)

        else:
            sigma = float(error['symerror'].strip('%')) / 100.0 * cv
            return sigma / np.sqrt(2.0), -sigma / np.sqrt(2.0)


def HEP_table_to_df_ATLAS_2JET_7TEV_R06(heptable, scenario='nominal'):
    """
    Given hep data table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and
    np.nan entries

    Parameters
    ----------
    heptable : str
            path to hepdata table

    scenario : 0, 1, 2
            0 = nominal, 1 = stronger, 2 = weaker
    """

    with open(heptable) as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][SCENARIO_ATLAS_2JET_7TEV_R06[scenario]]['values']
    df = pd.DataFrame(index=range(1, len(card) + 1))

    columns = {}
    errors = card[0]['errors']
    for i, err in enumerate(errors):
        # luminosity uncertainty, always symmetric
        if (
            (scenario == 'nominal' and i == 67)
            or (scenario == 'stronger' and i == 57)
            or (scenario == 'weaker' and i == 69)
        ):
            columns["lum"] = np.nan

        elif err['label'] == 'sys':
            columns[f"{err['label']}_plus_{i}"] = np.nan
            columns[f"{err['label']}_minus_{i}"] = np.nan

    df = pd.concat([df, pd.DataFrame(columns, index=df.index)], axis=1)
    return df


def fill_df_ATLAS_2JET_7TEV_R06(heptable, scenario='nominal'):
    """
    Fill a data frame with index
    corresponding to dijet mass bins
    and columns to different uncertainties
    Each df is for a fixed rapidity bin.
    """

    with open(heptable) as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][SCENARIO_ATLAS_2JET_7TEV_R06[scenario]]['values']
    df_nan = HEP_table_to_df_ATLAS_2JET_7TEV_R06(heptable, scenario)

    for i, dat in enumerate(card):
        cv = dat['value']

        for j, err in enumerate(dat['errors']):
            if (
                (scenario == 'nominal' and j == 67)
                or (scenario == 'stronger' and j == 57)
                or (scenario == 'weaker' and j == 69)
            ):
                # d_p, d_m = process_err_ATLAS_2JET_7TEV_R06(err,cv)
                # df_nan.loc[df_nan.index == i+1,"lum_plus"] = d_p
                # df_nan.loc[df_nan.index == i+1,"lum_minus"] = d_m
                err["lum"] = "lum"
                sigma = process_err_ATLAS_2JET_7TEV_R06(err, cv)

                df_nan.loc[df_nan.index == i + 1, "lum"] = sigma

            elif err['label'] == 'sys':
                d_p, d_m = process_err_ATLAS_2JET_7TEV_R06(err, cv)
                df_nan.loc[df_nan.index == i + 1, f"{err['label']}_plus_{j}"] = d_p
                df_nan.loc[df_nan.index == i + 1, f"{err['label']}_minus_{j}"] = d_m

    return df_nan
