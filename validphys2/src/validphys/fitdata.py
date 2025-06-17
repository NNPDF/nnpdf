"""
Utilities for loading data from fit folders
"""

from collections import OrderedDict, defaultdict, namedtuple
from dataclasses import dataclass
from io import StringIO
import json
import logging
import pathlib

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.checks import CheckError, make_argcheck
from reportengine.floatformatting import ValueErrorTuple
from reportengine.table import table
from validphys import checks
from validphys.core import PDF
from validphys.plotoptions.core import get_info
from validphys.utils import yaml_safe

# TODO: Add more stuff here as needed for postfit
LITERAL_FILES = ['chi2exps.log']
REPLICA_FILES = ['.dat', '.json']
FIT_SUMRULES = ["momentum", "uvalence", "dvalence", "svalence"]

# t = blessings.Terminal()
log = logging.getLogger(__name__)


def num_fitted_replicas(fit):
    """Function to obtain the number of nnfit replicas. That is
    the number of replicas before postfit was run.
    """
    with open(fit.path / "postfit" / "veto_count.json") as stream:
        veto = json.load(stream)
    # In principle we could use any keys
    return len(veto["Total"])


# TODO setup make_check on these
def check_nnfit_results_path(path):
    """Returns True if the requested path is a valid results directory,
    i.e if it is a directory and has a 'nnfit' subdirectory"""
    if not path.is_dir():
        log.warning(f"Path is not a directory {path}")
        return False
    if not (path / 'nnfit').is_dir():
        log.warning("Path {path/'nnfit'} is not a directory")
        return False
    return True


def check_lhapdf_info(results_dir, fitname):
    """Check that an LHAPDF info metadata file is
    present in the fit results"""
    info_path = results_dir.joinpath('nnfit', f'{fitname}.info')
    if not info_path.is_file():
        log.warning(f"Cannot find info file at {info_path}")
        return False
    return True


# TODO This should establish if the .dat files are corrupted or not
def check_replica_files(replica_path, prefix):
    """Verification of a replica results directory at `replica_path`
    for a fit named `prefix`. Returns True if the results
    directory is complete"""

    path = pathlib.Path(replica_path)
    if not path.is_dir():
        log.warning(f"Invalid directory for replica {path}")
        return False
    valid = True
    for f in LITERAL_FILES:
        test_path = path / f
        if not test_path.is_file():
            log.warning(f"Missing file: {test_path}")
            valid = False
    main_path = path / prefix
    for f in REPLICA_FILES:
        test_path = main_path.with_suffix(f)
        if not test_path.is_file():
            log.warning(f"Missing file: {test_path}")
            valid = False
    if not valid:
        log.warning(f"Found invalid replica {path}")
    return valid


@dataclass
class FitInfo:
    """Hold various metadata about the replicas of a fit

    Parameter
    ---------
    nite: best (stop) epoch of the replica
    training: Value of the training error function
    validation: Value of the validation error function
    chi2: Value of the experimental chi2
    pos_flag: State of the positivity pass flag
    arclengths: Array of the arc-length of each parton
    integrability: Array of values for the integrability checks
    """

    nite: int
    training: float
    validation: float
    chi2: float
    pos_flag: bool
    arclengths: np.ndarray
    integnumbers: np.ndarray

    @property
    def has_converged(self):
        """Uses the positivity flag as a proxy for convergence. Where convergence is defined mainly through the two constraints: the validation loss being smaller than the threshold and the positivity criterion being satisfied.
        When a fit does not reach convergence, the positivity flag is never flipped to true."""
        return self.pos_flag


def _old_load_fitinfo(old_fitinfo):
    """Process the data in the old ``.fitinfo`` files
    so that comparisons can still be run against very old fits
    """
    with old_fitinfo.open("r", encoding="utf-8") as fitinfo_file:
        fitinfo_line = fitinfo_file.readline().split()  # General fit properties
        fitinfo_arcl = fitinfo_file.readline()  # Replica arc-lengths
        fitinfo_integ = fitinfo_file.readline()  # Replica integ-numbers

        n_iterations = int(fitinfo_line[0])
        erf_validation = float(fitinfo_line[1])
        erf_training = float(fitinfo_line[2])
        chisquared = float(fitinfo_line[3])
        is_positive = fitinfo_line[4] == "POS_PASS"
        arclengths = np.fromstring(fitinfo_arcl, sep=' ')
        integnumbers = np.fromstring(fitinfo_integ, sep=' ')
    return FitInfo(
        n_iterations,
        erf_training,
        erf_validation,
        chisquared,
        is_positive,
        arclengths,
        integnumbers,
    )


def load_fitinfo(replica_path, prefix):
    """Process the data in the ``.json.`` file for a single replica into a ``FitInfo`` object.
    If the ``.json`` file does not exist an old-format fit is assumed and ``old_load_fitinfo``
    will be called instead.
    """
    p = replica_path / (prefix + ".json")
    if not p.exists():
        return _old_load_fitinfo(p.with_suffix(".fitinfo"))
    fitinfo_dict = json.loads(p.read_text(encoding="utf-8"))

    n_iterations = fitinfo_dict["best_epoch"]
    erf_validation = fitinfo_dict["erf_vl"]
    erf_training = fitinfo_dict["erf_tr"]
    chisquared = fitinfo_dict["chi2"]
    is_positive = fitinfo_dict["pos_state"] == "POS_PASS"
    arclengths = np.array(fitinfo_dict["arc_lengths"])
    integnumbers = np.array(fitinfo_dict["integrability"])
    return FitInfo(
        n_iterations,
        erf_training,
        erf_validation,
        chisquared,
        is_positive,
        arclengths,
        integnumbers,
    )


@checks.check_has_fitted_replicas
def replica_paths(fit):
    """Return the paths of all the replicas"""
    # Total number of members = number of replicas + 1
    l = len(PDF(fit.name))
    postfit_path = fit.path / 'postfit'
    old_postfit_path = fit.path / 'nnfit'
    if postfit_path.is_dir():
        return [postfit_path / f'replica_{index}' for index in range(1, l)]
    return [old_postfit_path / f'replica_{index}' for index in range(1, l)]


def replica_data(fit, replica_paths):
    """Load the necessary data from the ``.json`` file of each of the replicas.
    The corresponding PDF set must be installed in the LHAPDF path.

    The included information is:

    ('nite', 'training', 'validation', 'chi2', 'pos_status', 'arclenghts')"""
    return [load_fitinfo(path, fit.name) for path in replica_paths]


@table
def fit_summary(fit_name_with_covmat_label, replica_data, total_chi2_data, total_phi_data):
    """Summary table of fit properties
    - Central chi-squared
    - Average chi-squared
    - Training and Validation error functions
    - Training lengths
    - Phi

    Note:
    Chi-squared values from the replica_data are not used here (presumably
    they are fixed to being t0)

    This uses a corrected form for the error on phi in comparison to the
    vp1 value. The error is propagated from the uncertainty on the
    average chi-squared only.

    """
    nrep = len(replica_data)
    ndata = total_chi2_data.ndata
    central_chi2 = total_chi2_data.central_result / ndata
    member_chi2 = total_chi2_data.replica_result.error_members() / ndata

    nite = [x.nite for x in replica_data]
    etrain = [x.training for x in replica_data]
    evalid = [x.validation for x in replica_data]

    phi, _ = total_phi_data
    phi_err = np.std(member_chi2) / (2.0 * phi * np.sqrt(nrep))

    VET = ValueErrorTuple
    data = OrderedDict(
        (
            (r"$\chi^2$", f"{central_chi2:.5f}"),
            (r"$<E_{\mathrm{trn}}>$", f"{VET(np.mean(etrain), np.std(etrain))}"),
            (r"$<E_{\mathrm{val}}>$", f"{VET(np.mean(evalid), np.std(evalid))}"),
            (r"$<TL>$", f"{VET(np.mean(nite), np.std(nite))}"),
            (r"$<\chi^2>$", f"{VET(np.mean(member_chi2), np.std(member_chi2))}"),
            (r"$\phi$", f"{VET(phi, phi_err)}"),
        )
    )

    return pd.Series(data, index=data.keys(), name=fit_name_with_covmat_label)


collected_fit_summaries = collect('fit_summary', ('fits', 'fitcontext'))


@table
def summarise_fits(collected_fit_summaries):
    """Produces a table of basic comparisons between fits, includes
    all the fields used in fit_summary"""
    return pd.concat(collected_fit_summaries, axis=1)


@checks.check_use_t0
@table
def t0_chi2_info_table(pdf, dataset_inputs_abs_chi2_data, t0pdfset, use_t0):
    """Provides table with
    - t0pdfset name
    - Central t0-chi-squared
    - Average t0-chi-squared
    """
    ndata = dataset_inputs_abs_chi2_data.ndata
    central_chi2 = dataset_inputs_abs_chi2_data.central_result / ndata
    member_chi2 = dataset_inputs_abs_chi2_data.replica_result.error_members() / ndata

    VET = ValueErrorTuple
    data = OrderedDict(
        (
            ("t0pdfset", f"{t0pdfset}"),
            (r"$\chi^2_{t0}$", f"{central_chi2:.5f}"),
            (r"$<\chi^2_{t0}>$", f"{VET(np.mean(member_chi2), np.std(member_chi2))}"),
        )
    )

    return pd.Series(data, index=data.keys(), name=pdf.label)


fits_replica_data = collect('replica_data', ('fits',))

# Do collect in two parts so we get a list for each fit instead of a single list
all_datasets = collect('dataset', ('data',))
fits_datasets = collect('all_datasets', ('fits', 'fitinputcontext'))


@make_argcheck
def _assert_two_fits(fits):
    """Check that there are exatly two fits"""
    if len(fits) != 2:
        raise CheckError("Exactly two fits are required")


DatasetComp = namedtuple('DatasetComp', ('common', 'first_only', 'second_only'))


@_assert_two_fits
def match_datasets_by_name(fits, fits_datasets):
    """Return a tuple with common, first_only and second_only.
    The elements of the tuple are mappings where the keys are dataset names
    and the values are the two datasets contained in each fit for common, and
    the corresponfing dataset inclucded only in the first fit and only in the
    second fit."""

    firstds, secondds = ({ds.name: ds for ds in datasets} for datasets in fits_datasets)
    common_keys = firstds.keys() & secondds.keys()
    first_keys = firstds.keys() - secondds.keys()
    seccond_keys = secondds.keys() - firstds.keys()

    common = {k: (firstds[k], secondds[k]) for k in common_keys}
    first_only = {k: firstds[k] for k in first_keys}
    second_only = {k: secondds[k] for k in seccond_keys}
    return DatasetComp(common, first_only, second_only)


def _prepare_string(datalist, sort: bool = True):
    """Given a list of DataSetSpec, return a string with
    dataset_label (commondata.name)
    if sorted is true, it will be sorted by the name of the commondata"""
    if sort:
        datalist = sorted(datalist, key=str)
    str_list = []
    for dataset in datalist:
        info = get_info(dataset.commondata)
        str_list.append(f" - {info.dataset_label} (`{dataset}`)")
    return "\n".join(str_list)


def print_dataset_differences(fits, match_datasets_by_name, print_common: bool = True):
    """Given exactly two fits, print the datasets that are included in one "
    "but not in the other. If `print_common` is True, also print the datasets
    that are common.

    For the purposes of visual aid, everything is ordered by the dataset name,
    in terms of the the convention for the commondata means that everything is order by:
        1. Experiment name
        2. Process
        3. Energy
    """
    m = match_datasets_by_name
    first, second = fits
    res = StringIO()
    # When looking at the differences between dataset it is _very important_ to know also the theories
    th1 = first.as_input()["theory"]["theoryid"]
    th2 = second.as_input()["theory"]["theoryid"]
    if m.common and print_common:
        res.write(f"The following datasets are included in both `{first}` and `{second}`:\n\n")
        # Only one needed since they are common
        preprocess = [i[0] for i in m.common.values()]
        res.write(_prepare_string(preprocess))
        res.write("\n\n")
    if m.first_only:
        res.write(f"The following datasets are included in `{first}` but not in `{second}`:\n\n")
        res.write(_prepare_string(m.first_only.values()))
        res.write("\n\n")
    if m.second_only:
        res.write(f"The following datasets are included in `{second}` but not in `{first}`:\n\n")
        res.write(_prepare_string(m.second_only.values()))
        res.write("\n\n")

    if th1 != th2:
        res.write(f"`{first}` uses theoryid=`{th1}` while `{second}` used theoryid=`{th2}`\n\n")
    else:
        res.write(f"The theories used in the fits are identical (theoryid={th1})\n\n")
    if not first and not second and not print_common:
        res.write("The datasets included in the fits are identical.")

    return res.getvalue()


print_dataset_differences.highlight = 'markdown'


@_assert_two_fits
def test_for_same_cuts(fits, match_datasets_by_name):
    """Given two fits, return a list of tuples `(first, second)`
    where `first` and `second` are
    DatasetSpecs that correspond to the same dataset but have different cuts,
    such that `first` is included in the first fit and `second` in the second.
    """
    common = match_datasets_by_name.common
    first_fit, second_fit = fits
    res = []
    for ds, (first, second) in common.items():
        if first.cuts:
            c1 = first.cuts.load()
        else:
            c1 = np.arange(first.commondata.ndata)
        if second.cuts:
            c2 = second.cuts.load()
        else:
            c2 = np.arange(second.commondata.ndata)
        if not np.array_equal(c1, c2):
            msg = "Cuts for {} are not the same:\n{}:\n{}\n\n{}:\n{}".format(
                ds, first_fit, c1, second_fit, c2
            )
            log.info(msg)
            res.append((first, second))
    return res


def print_different_cuts(fits, test_for_same_cuts):
    """Print a summary of the datasets that are included in both fits but have
    different cuts."""
    res = StringIO()
    first_fit, second_fit = fits
    if test_for_same_cuts:
        res.write(
            "The following datasets are both included but have different kinematical cuts:\n\n"
        )
        for first, second in test_for_same_cuts:
            info = get_info(first.commondata)
            total_points = first.commondata.ndata
            res.write(" - %s:\n" % info.dataset_label)
            first_len = len(first.cuts.load()) if first.cuts else total_points
            second_len = len(second.cuts.load()) if second.cuts else total_points
            res.write(
                "    * %s includes %d out of %d points.\n" % (first_fit, first_len, total_points)
            )
            res.write(
                "    * %s includes %d out of %d points.\n" % (second_fit, second_len, total_points)
            )
        res.write('\n')

    return res.getvalue()


def fit_theory_covmat_summary(fit, fitthcovmat):
    """returns a table with a single column for the `fit`, with three rows
    indicating if the theory covariance matrix was used in the 'sampling' of the pseudodata,
    the 'fitting', and the 'validphys statistical estimators' in the current namespace for that fit.
    """
    try:
        config = fit.as_input()['theorycovmatconfig']
    except KeyError:
        config = {'use_thcovmat_in_sampling': False, 'use_thcovmat_in_fitting': False}
    sampling = config.get('use_thcovmat_in_sampling', False)
    fitting = config.get('use_thcovmat_in_fitting', False)
    report = bool(fitthcovmat)
    df = pd.DataFrame(
        [sampling, fitting, report],
        columns=[fit.label],
        index=['sampling', 'fitting', 'validphys statistical estimators'],
    )
    return df


fits_theory_covmat_summary = collect('fit_theory_covmat_summary', ('fits',))


@table
def summarise_theory_covmat_fits(fits_theory_covmat_summary):
    """Collects the theory covmat summary for all fits and concatenates them into a single table"""
    return pd.concat(fits_theory_covmat_summary, axis=1)


def _get_fitted_index(pdf, i):
    """Return the nnfit index for the replica i"""
    p = pdf.infopath.with_name(f'{pdf.name}_{i:04d}.dat')
    with open(p) as f:
        it = yaml_safe.load_all(f)
        metadata = next(it)
    return metadata['FromMCReplica']


@make_argcheck
def _check_has_replica_tags(pdf):
    """Check that the PDF has fitted index tags."""
    try:
        _get_fitted_index(pdf, 1)
    except KeyError as e:
        raise CheckError("PDF replica file don't contain the fitted replica tag.") from e


@_check_has_replica_tags
def fitted_replica_indexes(pdf):
    """Return nnfit index of replicas 1 to N."""
    return [_get_fitted_index(pdf, i) for i in range(1, len(pdf))]


fits_replica_indexes = collect('fitted_replica_indexes', ('fits', 'fitpdf'))


def fits_replica_data_correlated(fits_replica_data, fits_replica_indexes, fits):
    """Return a table with the same columns as ``replica_data`` indexed by the
    replica fit ID. For identical fits,
    the values across rows should be the same.

    If some replica ID is not present for a given fit (e.g. discarded by
    postfit), the corresponding entries in the table will be null.

    """
    dfs = []
    for dt, inds in zip(fits_replica_data, fits_replica_indexes):
        dfs.append(pd.DataFrame(dt, columns=FitInfo._fields, index=inds))
    return pd.concat(dfs, axis=1, keys=[fit.name for fit in fits])


@table
def datasets_properties_table(data_input):
    """Return dataset properties for each dataset in ``data_input``"""
    dataset_property_dict = defaultdict(list)
    for dataset in data_input:
        # only add elements if they don't evaluate to false
        ds_input_dict = {k: v for (k, v) in zip(dataset.argnames(), dataset.comp_tuple) if v}
        dataset_property_dict["Dataset"].append(ds_input_dict.pop("name"))
        dataset_property_dict["Training fraction"].append(ds_input_dict.pop("frac", "-"))
        dataset_property_dict["Weight"].append(ds_input_dict.pop("weight", "-"))
        dataset_property_dict["C-factors"].append(", ".join(ds_input_dict.pop("cfac", "-")))
        dataset_property_dict["Other fields"].append(
            ", ".join([f"{k}: {v}" for k, v in ds_input_dict.items()]) if ds_input_dict else "-"
        )
    df = pd.DataFrame(dataset_property_dict)
    df.set_index("Dataset", inplace=True)
    df = df[["Training fraction", "Weight", "C-factors", "Other fields"]]
    return df


@table
def fit_datasets_properties_table(fitinputcontext):
    """Returns table of dataset properties for each dataset used in a fit."""
    return datasets_properties_table(fitinputcontext["data_input"])


dataset_inputs_commondata = collect("commondata", ("data_input",))
groups_commondata = collect("dataset_inputs_commondata", ("group_dataset_inputs_by_metadata",))


def print_systype_overlap(groups_commondata, group_dataset_inputs_by_metadata):
    """Returns a set of systypes that overlap between groups.
    Discards the set of systypes which overlap but do not imply
    correlations

    """
    allow_list = {"CORR", "UNCORR", "SKIP", "THEORYUNCORR", "THEORYCORR"}
    systype_groups = dict()
    for group_cd, group in zip(groups_commondata, group_dataset_inputs_by_metadata):
        systype_groups[group["group_name"]] = {
            cd.load().systype_table.iloc[i]["name"]
            for cd in group_cd
            for i in range(cd.nsys)
            if cd.load().systype_table.iloc[i]["name"] not in allow_list
        }

    systype_overlap = set()
    groups_overlap = set()
    name_sys_list = list(systype_groups.items())
    for i, (name_i, sys_i) in enumerate(name_sys_list):
        for name_j, sys_j in name_sys_list[i + 1 :]:
            comp = sys_i.intersection(sys_j)
            if comp:
                groups_overlap.update({name_i, name_j})
            systype_overlap.update(comp)
    if systype_overlap:
        return groups_overlap, systype_overlap
    else:
        return "No overlap of systypes"


@table
def fit_code_version(fit):
    """Returns table with the code version from ``replica_1/{fitname}.json`` files.
    Note that the version for thensorflow distinguishes between the mkl=on and off version
    """
    vinfo = {}
    for json_path in fit.path.glob(f"nnfit/replica_*/{fit.name}.json"):
        tmp = json.loads(json_path.read_text(encoding="utf-8")).get("version")
        # if there's at least a replica without information, then they are all inconsistent
        if tmp is None:
            vinfo = {i: "inconsistent" for i in tmp}
            break
        elif not vinfo:
            vinfo = tmp
        for k, v in tmp.items():
            # If any value has changed for any replica, set it as inconsistent
            vinfo[k] = v if v == vinfo[k] else "inconsistent"

    return pd.DataFrame(vinfo.items(), columns=["module", fit.name]).set_index("module")


fits_fit_code_version = collect("fit_code_version", ("fits",))


@table
def fits_version_table(fits_fit_code_version):
    """Produces a table of version information for multiple fits."""
    vtable = pd.concat(fits_fit_code_version, axis=1)
    # Fill NaNs with "unavailable"
    vtable.fillna("unavailable", inplace=True)
    return vtable
