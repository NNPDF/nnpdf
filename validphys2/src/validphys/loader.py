"""
Resolve paths to useful objects, and query the existence of different resources
within the specified paths.
"""

import functools
import logging
import mimetypes
import os
import pathlib
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.parse as urls

import requests

from nnpdf_data import THEORY_CARDS_PATH
from nnpdf_data.commondataparser import parse_new_metadata, parse_set_metadata
from nnpdf_data.coredata import generate_path_filtered_data
from nnpdf_data.utils import get_nnpdf_profile
from nnpdf_data.validphys_compatibility import legacy_to_new_map, legacy_to_new_mapping
from reportengine import filefinder
from validphys import lhaindex
from validphys.core import (
    PDF,
    CommonDataSpec,
    Cuts,
    CutsPolicy,
    DataGroupSpec,
    DataSetSpec,
    FitSpec,
    FKTableSpec,
    HyperscanSpec,
    IntegrabilitySetSpec,
    InternalCutsWrapper,
    PositivitySetSpec,
    TheoryIDSpec,
)
from validphys.utils import tempfile_cleaner, yaml_safe

log = logging.getLogger(__name__)


class LoaderError(Exception):
    pass


class LoadFailedError(FileNotFoundError, LoaderError):
    pass


class DataNotFoundError(LoadFailedError):
    pass


class SysNotFoundError(LoadFailedError):
    pass


class FKTableNotFound(LoadFailedError):
    pass


class CfactorNotFound(LoadFailedError):
    pass


class CompoundNotFound(LoadFailedError):
    pass


class TheoryNotFound(LoadFailedError):
    pass


class EkoNotFound(LoadFailedError):
    pass


class TheoryMetadataNotFound(LoadFailedError):
    pass


class TheoryDataBaseNotFound(LoadFailedError):
    pass


class FitNotFound(LoadFailedError):
    pass


class HyperscanNotFound(LoadFailedError):
    pass


class CutsNotFound(LoadFailedError):
    pass


class PDFNotFound(LoadFailedError):
    pass


class ProfileNotFound(LoadFailedError):
    pass


class RemoteLoaderError(LoaderError):
    pass


class InconsistentMetaDataError(LoaderError):
    pass


def _fail_nicely_DataNotFoundError(setname):
    """Fail with a DataNotFoundError, but try to check whether the dataset exist in the
    translation layer, and if it does, offer the translation to the user."""
    new_name, _ = legacy_to_new_map(setname, None)
    err = ""
    if new_name != setname:
        err = f"\nNote that old names are no longer accepted. Perhaps you meant {new_name}"
    raise DataNotFoundError(f"Dataset {setname} not found. Is the name correct? {err}")


class LoaderBase:
    """
    Base class for the NNPDF loader.
    It can take as input a profile dictionary from which all data can be read.
    It is possible to override the datapath and resultpath when the class is instantiated.
    """

    def __init__(self, profile=None):
        if not isinstance(profile, dict):
            # If profile is a path, a str or None, read it from the default path
            profile = get_nnpdf_profile(profile)

        # Retrieve important paths from the profile if not given
        datapaths = [pathlib.Path(i) for i in profile["data_path"]]
        theories_path = pathlib.Path(profile["theories_path"])
        resultspath = pathlib.Path(profile["results_path"])
        ekos_path = pathlib.Path(profile["ekos_path"])

        # Create the theories and results paths if they don't exist already
        theories_path.mkdir(exist_ok=True, parents=True)
        ekos_path.mkdir(exist_ok=True, parents=True)
        resultspath.mkdir(exist_ok=True, parents=True)

        # And save them up
        self.commondata_folders = tuple(datapaths)
        self._theories_path = theories_path
        self._ekos_path = ekos_path
        self.resultspath = resultspath
        self._extremely_old_fits = set()
        self.nnprofile = profile

    @property
    def hyperscan_resultpath(self):
        hyperscan_path = pathlib.Path(self.nnprofile["hyperscan_path"])
        hyperscan_path.mkdir(parents=True, exist_ok=True)
        return hyperscan_path

    def _vp_cache(self):
        """Return the vp-cache path, and create it if it doesn't exist"""
        vpcache = pathlib.Path(self.nnprofile['validphys_cache_path'])
        if not vpcache.exists():
            try:
                log.info(f"Creating validphys cache directory: {vpcache}")
                vpcache.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise LoaderError("Could not create the cache directory " f"at {vpcache}") from e
        return vpcache


class Loader(LoaderBase):
    """Load various resources from the NNPDF data path."""

    @property
    def available_fits(self):
        try:
            return [p.name for p in self.resultspath.iterdir() if p.is_dir()]
        except OSError:
            return []

    @property
    def available_hyperscans(self):
        try:
            return [p.name for p in self.hyperscan_resultpath.iterdir() if p.is_dir()]
        except OSError:
            return []

    @property
    @functools.lru_cache
    def available_theories(self):
        """Return a string token for each of the available theories"""
        theory_token = 'theory_'
        return {
            folder.name[len(theory_token) :]
            for folder in self._theories_path.glob(theory_token + '*')
        }

    @property
    @functools.lru_cache
    def available_ekos(self):
        """Return a string token for each of the available theories"""
        return {
            eko_path.parent.name.split("_")[1] for eko_path in self._theories_path.glob("*/eko.tar")
        }

    @property
    @functools.lru_cache
    def available_datasets(self):
        """Provide all available datasets that were available before the new commondata
        was implemented and that have a translation.
        Returns old names
        """
        # TODO: return new names... without loading the metadata...
        # Skip Positivity and Integrability
        skip = ["POS", "INTEG"]
        # Skip datasets for which a translation exists but were not implemented in the old way
        skip += ["STAR", "ATLAS_WJ_JET_8TEV_"]
        # skip all the _DW_ datasets that are repeated
        skip += [i for i in legacy_to_new_mapping.keys() if "_DW_" in i]
        old_datasets = [i for i in legacy_to_new_mapping.keys() if not i.startswith(tuple(skip))]
        return set(old_datasets)

    @property
    @functools.lru_cache
    def implemented_datasets(self):
        """Provide all implemented datasets that can be found in the datafiles folder
        regardless of whether they can be used for fits (i.e., whether they include a theory),
        are "fake" (integrability/positivity) or are missing some information.
        """
        datasets = []
        for commondata_folder in self.commondata_folders:
            for metadata_file in commondata_folder.glob("*/metadata.yaml"):
                datasets += parse_set_metadata(metadata_file).allowed_datasets
        return datasets

    @property
    @functools.lru_cache
    def available_pdfs(self):
        return lhaindex.expand_local_names('*')

    def check_commondata(
        self, setname, sysnum=None, use_fitcommondata=False, fit=None, variant=None
    ):
        """Prepare the commondata files to be loaded.
        A commondata is defined by its name (``setname``) and the variant(s) (``variant``)

        The function ``parse_dataset_input`` in ``config.py`` translates all known old commondata
        into their new names (and variants),
        therefore this function should only receive requestes for new format.

        Any actions trying to requests an old-format commondata from this function will log
        an error message. This error message will eventually become an actual error.
        """
        if use_fitcommondata:
            if not fit:
                raise LoadFailedError("Must specify a fit when setting use_fitcommondata")

            # First, load the base commondata which will be used as container
            basedata = self.check_commondata(setname, variant=variant, sysnum=sysnum)
            # and the filename for the new data
            data_path, unc_path = generate_path_filtered_data(fit.path, setname)

            if not data_path.exists():
                # We might be dealing with legacy names and with legacy paths
                return DataNotFoundError(f"No fit data found for {setname} ({data_path})")

            return basedata.with_modified_data(data_path, uncertainties_file=unc_path)

        # Get data folder and observable name and check for existence
        try:
            setfolder, observable_name = setname.rsplit("_", 1)
        except ValueError:
            _fail_nicely_DataNotFoundError(setname)

        for commondata_folder in self.commondata_folders:
            set_path = commondata_folder / setfolder
            if set_path.exists():
                break
        else:
            _fail_nicely_DataNotFoundError(setname)

        metadata_path = set_path / "metadata.yaml"
        metadata = parse_new_metadata(metadata_path, observable_name, variant=variant)
        return CommonDataSpec(setname, metadata)

    @functools.lru_cache
    def check_theoryID(self, theoryID):
        theoryID = int(theoryID)
        theopath = self._theories_path / f"theory_{theoryID}"
        if not theopath.exists():
            raise TheoryNotFound(f"Could not find theory {theoryID}. Folder '{theopath}' not found")
        return TheoryIDSpec(theoryID, theopath, self.theorydb_folder)

    @functools.lru_cache
    def check_eko(self, theoryID):
        """Check the eko exists and return the path to it"""
        eko_path = self._ekos_path / f"eko_{int(theoryID)}.tar"
        if not eko_path.exists():
            raise EkoNotFound(f"Could not find eko {eko_path} in theory: {theoryID}")
        return eko_path

    @property
    def theorydb_folder(self):
        """Checks theory db file exists and returns path to it"""
        if THEORY_CARDS_PATH.exists():
            return THEORY_CARDS_PATH
        raise TheoryDataBaseNotFound(
            f"could not find theory db folder. Directory not found at {THEORY_CARDS_PATH}"
        )

    def check_fktable(self, theoryID, setname, cfac):
        _, theopath = self.check_theoryID(theoryID)
        fkpath = theopath / 'fastkernel' / ('FK_%s.dat' % setname)
        if not fkpath.exists():
            raise FKTableNotFound(
                f"Could not find FKTable for set '{setname}'. File '{fkpath}' not found"
            )

        cfactors = self.check_cfactor(theoryID, setname, cfac)
        return FKTableSpec(fkpath, cfactors)

    def check_fk_from_theory_metadata(self, theory_metadata, theoryID, cfac=None):
        """Load a pineappl fktable in the new commondata forma
        Receives a theory metadata describing the fktables necessary for a given observable
        the theory ID and the corresponding cfactors.
        The cfactors should correspond directly to the fktables, the "compound folder"
        is not supported for pineappl theories. As such, the name of the cfactor is expected to be
            CF_{cfactor_name}_{fktable_name}
        """
        theory = self.check_theoryID(theoryID)
        fklist = theory_metadata.fktables_to_paths(theory.path / "fastkernel")
        op = theory_metadata.operation

        if not cfac or cfac is None:
            fkspecs = [FKTableSpec(i, None, theory_metadata) for i in fklist]
            return fkspecs, op

        cfactors = []
        for operand in theory_metadata.FK_tables:
            tmp = [self.check_cfactor(theoryID, fkname, cfac) for fkname in operand]
            cfactors.append(tuple(tmp))

        fkspecs = [FKTableSpec(i, c, theory_metadata) for i, c in zip(fklist, cfactors)]
        return fkspecs, theory_metadata.operation

    def check_compound(self, theoryID, setname, cfac):
        thid, theopath = self.check_theoryID(theoryID)
        compound_spec_path = theopath / 'compound' / ('FK_%s-COMPOUND.dat' % setname)
        try:
            with compound_spec_path.open() as f:
                # Drop first line with comment
                next(f)
                txt = f.read()
        except FileNotFoundError as e:
            msg = "Could not find COMPOUND set '%s' for theory %d: %s" % (setname, int(thid), e)
            raise CompoundNotFound(msg)
        # This is a little bit funny, but is the least amount of thinking...
        yaml_format = 'FK:\n' + re.sub('FK:', ' - ', txt)
        data = yaml_safe.load(yaml_format)
        # we have to split out 'FK_' the extension to get a name consistent
        # with everything else
        try:
            tables = [self.check_fktable(theoryID, name[3:-4], cfac) for name in data['FK']]
        except FKTableNotFound as e:
            raise LoadFailedError(
                f"Incorrect COMPOUND file '{compound_spec_path}'. "
                f"Searching for non-existing FKTable:\n{e}"
            ) from e
        op = data['OP']
        return tuple(tables), op

    def check_cfactor(self, theoryID, setname, cfactors):
        _, theopath = self.check_theoryID(theoryID)
        cf = []
        for cfactor in cfactors:
            cfactorpath = theopath / "cfactor" / f"CF_{cfactor}_{setname}.dat"
            if not cfactorpath.exists():
                msg = (
                    f"Could not find cfactor '{cfactor}' for FKTable {setname}."
                    f"File {cfactorpath} does not exist in {theoryID}"
                )
                raise CfactorNotFound(msg)
            cf.append(cfactorpath)

        return tuple(cf)

    def _check_lagrange_multiplier_set(self, theoryID, setname):
        """Check an integrability or positivity dataset"""
        cd = self.check_commondata(setname, 'DEFAULT')
        th = self.check_theoryID(theoryID)
        fk, _ = self._check_theory_old_or_new(th, cd, [])
        return cd, fk, th

    def check_posset(self, theoryID, setname, postlambda, rules):
        """Load a positivity dataset"""
        cd, fk, th = self._check_lagrange_multiplier_set(theoryID, setname)
        return PositivitySetSpec(setname, cd, fk, postlambda, th, rules)

    def check_integset(self, theoryID, setname, postlambda, rules):
        """Load an integrability dataset"""
        cd, fk, th = self._check_lagrange_multiplier_set(theoryID, setname)
        return IntegrabilitySetSpec(setname, cd, fk, postlambda, th, rules)

    def check_fit(self, fitname):
        resultspath = self.resultspath
        if fitname != pathlib.Path(fitname).name:
            raise FitNotFound(
                f"Could not find fit '{fitname}' in '{resultspath} "
                "because the name doesn't correspond to a valid filename"
            )
        p = resultspath / fitname
        if p.is_dir():
            return FitSpec(fitname, p)
        if not p.is_dir():
            msg = f"Could not find fit '{fitname}' in '{resultspath}'. Folder '{p}' not found"
            raise FitNotFound(msg)
        msg = f"Could not load fit '{fitname}' from '{resultspath}. '{p}' must be a folder"
        raise FitNotFound(msg)

    def check_hyperscan(self, hyperscan_name):
        """Obtain a hyperscan run"""
        resultspath = self.hyperscan_resultpath
        if hyperscan_name != pathlib.Path(hyperscan_name).name:
            raise HyperscanNotFound(
                f"Could not find fit '{hyperscan_name}' in '{resultspath} "
                "because the name doesn't correspond to a valid filename"
            )
        p = resultspath / hyperscan_name
        if p.is_dir():
            hyperspec = HyperscanSpec(hyperscan_name, p)
            if hyperspec.tries_files:
                return hyperspec
            raise HyperscanNotFound(f"No hyperscan output find in {hyperscan_name}")

        raise HyperscanNotFound(
            f"Could not find hyperscan '{hyperscan_name}' in '{resultspath}'."
            f" Folder '{hyperscan_name}' not found"
        )

    def check_default_filter_rules(self, theoryid, defaults=None):
        # avoid circular import
        from validphys.filters import (
            Rule,
            default_filter_rules_input,
            default_filter_settings_input,
        )

        th_params = theoryid.get_description()
        if defaults is None:
            defaults = default_filter_settings_input()
        return tuple(
            Rule(inp, defaults=defaults, theory_parameters=th_params, loader=self)
            for inp in default_filter_rules_input()
        )

    def _check_theory_old_or_new(self, theoryid, commondata, cfac):
        """Given a theory and a commondata and a theory load the right fktable
        checks whether:
            1. the theory is a pineappl theory
            2. Select the right information (commondata name, legacy name or theory meta)
        """
        theoryno, _ = theoryid
        name = commondata.name
        if theoryid.is_pineappl():
            if (thmeta := commondata.metadata.theory) is None:
                # Regardless of the type of theory, request the existence of the field
                raise TheoryMetadataNotFound(f"No theory metadata found for {name}")
            fkspec, op = self.check_fk_from_theory_metadata(thmeta, theoryno, cfac)
        elif commondata.legacy_names is None:
            raise LoadFailedError(f"Cannot use an old theory with a purely new dataset ({name})")
        else:
            # Old theories can only be used with datasets that have a corresponding
            # old name to map to, and so we need to be able to load the cd at this point
            for legacy_name in commondata.legacy_names:
                # This might be slow, if it becomes a problem, the map function can be used instead
                try:
                    fkspec, op = self.check_compound(theoryno, legacy_name, cfac)
                except CompoundNotFound:
                    fkspec = self.check_fktable(theoryno, legacy_name, cfac)
                    op = None
                except FKTableNotFound:
                    continue
                break
            else:
                # If no fktable has been found after looping over legacy names, raise an Error
                raise FKTableNotFound(
                    "Could not find old FKtable or CompoundFile for any of the legacy names of {name} {commondata.legacy_names}"
                )
        return fkspec, op

    def check_dataset(
        self,
        name,
        *,
        rules=None,
        sysnum=None,
        theoryid,
        cfac=(),
        frac=1,
        cuts=CutsPolicy.INTERNAL,
        use_fitcommondata=False,
        fit=None,
        weight=1,
        variant=None,
    ):
        """Loads a given dataset
        If the dataset contains new-type fktables, use the
        pineappl loading function, otherwise fallback to legacy
        """
        if not isinstance(theoryid, TheoryIDSpec):
            theoryid = self.check_theoryID(theoryid)

        # TODO:
        # The dataset is checked twice, once here
        # and once by config in produce_commondata
        # once of the two __must__ be superfluous
        # note that both use information from dataset_input
        commondata = self.check_commondata(
            name, sysnum, use_fitcommondata=use_fitcommondata, fit=fit, variant=variant
        )

        fkspec, op = self._check_theory_old_or_new(theoryid, commondata, cfac)

        # Note this is simply for convenience when scripting. The config will
        # construct the actual Cuts object by itself
        if isinstance(cuts, str):
            cuts = CutsPolicy(cuts)
        if isinstance(cuts, CutsPolicy):
            if cuts is CutsPolicy.NOCUTS:
                cuts = None
            elif cuts is CutsPolicy.FROMFIT:
                cuts = self.check_fit_cuts(commondata, fit)
            elif cuts is CutsPolicy.INTERNAL:
                if rules is None:
                    rules = self.check_default_filter_rules(theoryid)
                cuts = self.check_internal_cuts(commondata, rules)
            elif cuts is CutsPolicy.FROM_CUT_INTERSECTION_NAMESPACE:
                raise LoaderError(f"Intersection cuts not supported in loader calls.")

        return DataSetSpec(
            name=name,
            commondata=commondata,
            fkspecs=fkspec,
            thspec=theoryid,
            cuts=cuts,
            frac=frac,
            op=op,
            weight=weight,
            rules=rules,
        )

    def check_experiment(self, name: str, datasets: list[DataSetSpec]) -> DataGroupSpec:
        """Loader method for instantiating DataGroupSpec objects. The NNPDF::Experiment
        object can then be instantiated using the load method.

        Parameters
        ----------
        name: str
            A string denoting the name of the resulting DataGroupSpec object.
        dataset: List[DataSetSpec]
            A list of DataSetSpec objects pre-created by the user. Note, these too
            will be loaded by Loader.

        Returns
        -------
        DataGroupSpec

        Example
        -------
        >>> from validphys.loader import Loader
        >>> l = Loader()
        >>> ds = l.check_dataset("CDF_Z0_1P96TEV_ZRAP", theoryid=40_000_000, cuts="internal")
        >>> exp = l.check_experiment("My DataGroupSpec Name", [ds])
        """
        if not isinstance(datasets, list):
            raise TypeError("Must specify a list of DataSetSpec objects to use")

        return DataGroupSpec(name, datasets)

    def check_pdf(self, name):
        if lhaindex.isinstalled(name):
            return PDF(name)
        raise PDFNotFound(name)

    def check_fit_cuts(self, commondata, fit):
        setname = commondata.name
        if fit is None:
            raise TypeError("Must specify a fit to use the cuts.")
        if not isinstance(fit, FitSpec):
            fit = self.check_fit(fit)
        _, fitpath = fit

        cuts_path = (fitpath / 'filter') / setname / ('FKMASK_' + setname + '.dat')

        # After 4.0.9 we changed to a new commondata format
        # In order to utilize cuts from old fits in new fits it is necessary to translate the names
        # There are two translation that might be necessary:
        # 1. New names in the runcard, old cuts in the 'fromfit' fit
        # 2. Old names in the runcard, new cuts in the 'fromfit' fit
        # In order to enforce the usage of the new names, only (1.) is implemented

        if not cuts_path.parent.exists():
            if commondata.legacy_names is None:
                raise CutsNotFound(f"Bad filter configuration. Could not find {cuts_path.parent}")

            # Maybe we have an old-name version of the cuts?
            for old_name in commondata.legacy_names:
                # Then, check whether there are cuts with the corresponding old name
                old_dir = cuts_path.parent.with_name(old_name)
                if old_dir.exists():
                    cuts_path = old_dir / f"FKMASK_{old_name}.dat"
                    break
            else:
                raise CutsNotFound(f"Bad filter configuration. Could not find {cuts_path.parent}")

        if not cuts_path.exists():
            cuts_path = None
        return Cuts(commondata, cuts_path)

    def check_internal_cuts(self, commondata, rules):
        return InternalCutsWrapper(commondata, rules)

    def check_vp_output_file(self, filename, extra_paths=('.',)):
        """Find a file in the vp-cache folder, or (with higher priority) in
        the ``extra_paths``."""
        try:
            vpcache = self._vp_cache()
        except KeyError as e:
            log.warning("Entry validphys_cache_path expected but not found in the nnprofile.")
        else:
            extra_paths = (*extra_paths, vpcache)

        finder = filefinder.FallbackFinder(extra_paths)
        try:
            path, name = finder.find(filename)
        except FileNotFoundError as e:
            raise LoadFailedError(f"Could not find '{filename}'") from e
        except filefinder.FinderError as e:
            raise LoaderError(e) from e
        return path / name


# http://stackoverflow.com/a/15645088/1007990
def _download_and_show(response, stream):
    total_length = response.headers.get('content-length')

    if total_length is None or not log.isEnabledFor(logging.INFO):
        stream.write(response.content)
    else:
        dl = 0
        prev_done = -1
        total_length = int(total_length)
        for data in response.iter_content(chunk_size=4096):
            dl += len(data)
            stream.write(data)
            if sys.stdout.isatty():
                done = int(50 * dl / total_length)
                if prev_done != done:
                    sys.stdout.write(f"\r[{'=' * done}{' '*(50 - done)}] ({done * 2}%)")
                    prev_done = done
                    sys.stdout.flush()
        sys.stdout.write('\n')


def download_file(url, stream_or_path, make_parents=False, delete_on_failure=False):
    """Download a file and show a progress bar if the INFO log level is
    enabled. If ``make_parents`` is ``True`` ``stream_or_path``
    is path-like, all the parent folders will
    be created."""
    # There is a bug in CERN's
    # Apache that incorrectly sets the Content-Encodig header to gzip, even
    # though it doesn't compress two times.
    # See: http://mail-archives.apache.org/mod_mbox/httpd-dev/200207.mbox/%3C3D2D4E76.4010502@talex.com.pl%3E
    # and e.g. https://bugzilla.mozilla.org/show_bug.cgi?id=610679#c30
    # If it looks like the url is already encoded, we do not request
    # it to be compressed
    headers = {}
    if mimetypes.guess_type(url)[1] is not None:
        headers['Accept-Encoding'] = None

    response = requests.get(url, stream=True, headers=headers)

    response.raise_for_status()

    if isinstance(stream_or_path, (str, bytes, os.PathLike)):
        p = pathlib.Path(stream_or_path)
        if p.is_dir():
            raise IsADirectoryError(p)
        log.info("Downloading %s to %s.", url, stream_or_path)
        if make_parents:
            p.parent.mkdir(exist_ok=True, parents=True)

        with tempfile.NamedTemporaryFile(
            delete=delete_on_failure, dir=p.parent, prefix=p.name, suffix='.part'
        ) as f:
            _download_and_show(response, f)
            shutil.move(f.name, p)
    else:
        log.info("Downloading %s.", url)
        _download_and_show(response, stream_or_path)


def download_and_extract(url, local_path, target_name=None):
    """Download a compressed archive and then extract it to the given path"""
    local_path = pathlib.Path(local_path)
    if not local_path.is_dir():
        raise NotADirectoryError(local_path)
    name = url.split('/')[-1]
    archive_dest = tempfile.NamedTemporaryFile(delete=False, suffix=name, dir=local_path)
    with archive_dest as t:
        log.debug("Saving data to %s", t.name)
        download_file(url, t)
    log.info("Extracting archive to %s", local_path)
    try:
        with tarfile.open(archive_dest.name) as res_tar:
            # Extract to a temporary directory
            folder_dest = tempfile.TemporaryDirectory(dir=local_path, suffix=name)
            dest_path = pathlib.Path(folder_dest.name)
            try:
                res_tar.extractall(path=dest_path, filter="data")
            except TypeError as e:
                if sys.version_info > (3, 9, 16):
                    # Filter was added in 3.9.17, and raises TypeError before then
                    # mac's default is still in 3.9.6, so fallback to the unsafe behaviour
                    raise e
                res_tar.extractall(path=dest_path)
            except tarfile.LinkOutsideDestinationError as e:
                if sys.version_info > (3, 11):
                    raise e
                # For older versions of python ``filter=data`` might be too restrictive
                # for the links inside the ``postfit`` folder if you are using more than one disk
                res_tar.extractall(path=dest_path, filter="tar")

            # Check there are no more than one item in the top level
            top_level_stuff = list(dest_path.glob("*"))
            if len(top_level_stuff) > 1:
                raise RemoteLoaderError(f"More than one item in the top level directory of {url}")

            if target_name is None:
                target_path = local_path
            else:
                target_path = local_path / target_name
            shutil.move(top_level_stuff[0], target_path)

    except Exception as e:
        log.error(
            f"The original archive at {t.name} was only extracted partially at \n{local_path}"
        )
        raise e
    else:
        os.unlink(archive_dest.name)


def _key_or_loader_error(f):
    @functools.wraps(f)
    def f_(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except KeyError as e:
            log.error(f"nnprofile is configured improperly: Key {e} is missing from the profile!")
            raise LoaderError(
                "Cannot attempt download because "
                "nnprofile is configured improperly: "
                f"Missing key '{e}'"
            ) from e

    return f_


# TODO: Make this async someday
class RemoteLoader(LoaderBase):
    @property
    @_key_or_loader_error
    def fit_urls(self):
        return self.nnprofile['fit_urls']

    @property
    @_key_or_loader_error
    def fit_index(self):
        return self.nnprofile['fit_index']

    @property
    @_key_or_loader_error
    def hyperscan_url(self):
        return self.nnprofile['hyperscan_urls']

    @property
    @_key_or_loader_error
    def hyperscan_index(self):
        return self.nnprofile['hyperscan_index']

    @property
    @_key_or_loader_error
    def theory_urls(self):
        return self.nnprofile['theory_urls']

    @property
    @_key_or_loader_error
    def theory_index(self):
        return self.nnprofile['theory_index']

    @property
    @_key_or_loader_error
    def eko_index(self):
        return self.nnprofile['eko_index']

    @property
    @_key_or_loader_error
    def eko_urls(self):
        return self.nnprofile['eko_urls']

    @property
    @_key_or_loader_error
    def nnpdf_pdfs_urls(self):
        return self.nnprofile['nnpdf_pdfs_urls']

    @property
    @_key_or_loader_error
    def nnpdf_pdfs_index(self):
        return self.nnprofile['nnpdf_pdfs_index']

    @property
    @_key_or_loader_error
    def lhapdf_urls(self):
        urls = self.nnprofile['lhapdf_urls']
        if len(urls) > 1:
            log.warning("Only one lhapdf_url is supported at the moment.")
        if len(urls) == 0:
            raise LoaderError("The specification for lhapdf_urls is empty in nnprofile")
        return urls

    def _remote_files_from_url(self, url, index, thing='files'):
        index_url = url + index
        try:
            resp = requests.get(index_url)
            resp.raise_for_status()
        except Exception as e:
            raise RemoteLoaderError(f"Failed to fetch remote {thing} index {index_url}: {e}") from e

        try:
            info = resp.json()['files']
        except Exception as e:
            raise RemoteLoaderError(
                f"Malformed index {index_url}. Expecting json with a key 'files': {e}"
            ) from e

        return {file.split('.')[0]: url + file for file in info}

    def remote_files(self, urls, index, thing='files'):
        d = {}
        for url in urls:
            try:
                d.update(self._remote_files_from_url(url, index, thing))
            except RemoteLoaderError as e:
                log.error(e)
        return d

    @property
    @functools.lru_cache
    def remote_fits(self):
        return self.remote_files(self.fit_urls, self.fit_index, thing="fits")

    @property
    @functools.lru_cache
    def remote_hyperscans(self):
        return self.remote_files(self.hyperscan_url, self.hyperscan_index, thing="hyperscan")

    @property
    @functools.lru_cache
    def remote_theories(self):
        token = 'theory_'
        rt = self.remote_files(self.theory_urls, self.theory_index, thing="theories")
        return {k[len(token) :]: v for k, v in rt.items()}

    @property
    @functools.lru_cache
    def remote_ekos(self):
        token = 'eko_'
        rt = self.remote_files(self.eko_urls, self.eko_index, thing="ekos")
        return {k[len(token) :]: v for k, v in rt.items()}

    @property
    @functools.lru_cache
    def remote_nnpdf_pdfs(self):
        return self.remote_files(self.nnpdf_pdfs_urls, self.nnpdf_pdfs_index, thing="PDFs")

    @functools.cached_property
    def remote_keywords(self):
        root = self.nnprofile['reports_root_url']
        url = urls.urljoin(root, 'index.json')
        try:
            req = requests.get(url)
            req.raise_for_status()
            keyobjs = req.json()['keywords']
            l = [k[0] for k in keyobjs]
        except requests.RequestException as e:
            raise RemoteLoaderError(e) from e
        return l

    @property
    def downloadable_fits(self):
        return list(self.remote_fits)

    @property
    def downloadable_hyperscans(self):
        return list(self.remote_hyperscans)

    @property
    def downloadable_theories(self):
        return list(self.remote_theories)

    @property
    def downloadable_ekos(self):
        return list(self.remote_ekos)

    @property
    def lhapdf_pdfs(self):
        return lhaindex.expand_index_names('*')

    @property
    def nnpdf_pdfs(self):
        return list(self.remote_nnpdf_pdfs)

    @property
    def downloadable_pdfs(self):
        return {*self.lhapdf_pdfs, *self.downloadable_fits, *self.nnpdf_pdfs}

    def download_fit(self, fitname):
        if not fitname in self.remote_fits:
            raise FitNotFound(f"Could not find fit '{fitname}' in remote index {self.fit_index}")

        with tempfile_cleaner(
            root=self.resultspath,
            exit_func=shutil.rmtree,
            exc=KeyboardInterrupt,
            prefix='fit_download_deleteme_',
        ) as tempdir:
            download_and_extract(self.remote_fits[fitname], tempdir)
            # Handle old-style fits compressed with 'results' as root.
            old_style_res = tempdir / 'results'
            if old_style_res.is_dir():
                move_target = old_style_res / fitname
            else:
                move_target = tempdir / fitname
            if not move_target.is_dir():
                raise RemoteLoaderError(
                    f"Unknown format for fit in {tempdir}. Expecting a folder {move_target}"
                )

            fitpath = self.resultspath / fitname
            shutil.move(move_target, fitpath)

        if lhaindex.isinstalled(fitname):
            log.warning(
                f"The PDF corresponding to the downloaded fit '{fitname}' "
                "exists in the LHAPDF path."
                " Will be erased and replaced with the new one."
            )
            p = pathlib.Path(lhaindex.finddir(fitname))
            if p.is_symlink():
                p.unlink()
            else:
                shutil.rmtree(p)
        else:
            p = pathlib.Path(lhaindex.get_lha_datapath()) / fitname
            # This is needed here as well because the path may be a
            # broken symlink.
            if p.is_symlink():
                p.unlink()
        gridpath = fitpath / 'postfit' / fitname
        gridpath_old = fitpath / 'nnfit' / fitname
        if gridpath.is_dir():
            p.symlink_to(gridpath, target_is_directory=True)
        else:
            log.warning(f"Cannot find {gridpath}. Falling back to old behaviour")
            p.symlink_to(gridpath_old, target_is_directory=True)

    def download_hyperscan(self, hyperscan_name):
        """Download a hyperscan run from the remote server
        Downloads the run to the results folder
        """
        if not hyperscan_name in self.remote_hyperscans:
            raise HyperscanNotFound(
                f"Could not find hyperscan {hyperscan_name} in remote index {self.hyperscan_index}"
            )

        with tempfile_cleaner(
            root=self.hyperscan_resultpath,
            exit_func=shutil.rmtree,
            exc=KeyboardInterrupt,
            prefix='fit_download_deleteme_',
        ) as tempdir:
            download_and_extract(self.remote_hyperscans[hyperscan_name], tempdir)
            move_target = tempdir / hyperscan_name
            if not move_target.is_dir():
                raise RemoteLoaderError(
                    f"Unknown format for fit in {tempdir}. Expecting a folder {move_target}"
                )
            hyperscan_path = self.hyperscan_resultpath / hyperscan_name
            shutil.move(move_target, hyperscan_path)

    def download_pdf(self, name):
        # Check if the pdf is an existing fit first
        try:
            # We don't want to download the fit here
            fit = Loader.check_fit(self, name)
        except FitNotFound:
            pass
        else:
            p = pathlib.Path(lhaindex.get_lha_datapath()) / fit.name
            fitpath = fit.path / 'postfit'
            fitpath_old = fit.path / 'nnfit'
            if fitpath.exists() or fitpath_old.exists():
                log.info(
                    "Found existing fit with the same name as the "
                    "requested PDF (%s). Symlinking the grid to the LHAPDF path (%s).",
                    name,
                    p,
                )
                # This is needed here as well because the path may be a
                # broken symlink.
                if p.is_symlink():
                    p.unlink()
                if fitpath.exists():
                    p.symlink_to(fitpath / fit.name)
                else:
                    p.symlink_to(fitpath_old / fit.name)
                return

        # It would be good to use the LHAPDF command line, except that it does
        # questionable things like returning 0 exit status when it fails to
        # download.
        _saved_exception = False
        if name in self.lhapdf_pdfs:
            try:
                url = self.lhapdf_urls[0] + name + '.tar.gz'
                # url = 'https://data.nnpdf.science/thisisatesttodelete/NNPDF31_nlo_as_0118.tar.gz'
                # url = 'https://data.nnpdf.science/patata/NNPDF31_nlo_as_0118.tar.gz'
                return download_and_extract(url, lhaindex.get_lha_datapath())
            except shutil.ReadError as e:
                _saved_exception = e
                log.error(
                    f"{e}. It seems the LHAPDF URLs aren't behaving, "
                    f"attempting to find resource in other repositories"
                )
                pass
            except requests.RequestException as e:
                _saved_exception = e
                log.error(
                    f"There was a problem with the connection: {e}. "
                    f"Attempting to find resource elsewhere."
                )
                pass
            except RemoteLoaderError as e:
                _saved_exception = e
                log.error(f"Failed to download resource: {e}. Attempting " f"to find it elsewhere.")
                pass
        if name in self.downloadable_fits:
            try:
                return self.download_fit(name)
            except requests.RequestException as e:
                _saved_exception = e
                log.error(
                    f"There was a problem with the connection: {e}. "
                    f"Attempting to find resource elsewhere."
                )
                pass
            except RemoteLoaderError as e:
                _saved_exception = e
                log.error(f"Failed to download resource: {e}. Attempting " f"to find it elsewhere.")
                pass
        if name in self.remote_nnpdf_pdfs:
            return download_and_extract(self.remote_nnpdf_pdfs[name], lhaindex.get_lha_datapath())
        elif _saved_exception:
            raise LoadFailedError(
                f"{_saved_exception}. The resource could not " f"be found elsewhere."
            ) from _saved_exception
        else:
            raise PDFNotFound("PDF '%s' is neither an uploaded fit nor an " "LHAPDF set." % name)

    def download_theoryID(self, thid):
        thid = str(int(thid))
        remote = self.remote_theories
        if thid not in remote:
            raise TheoryNotFound("Theory %s not available." % thid)
        download_and_extract(remote[thid], self._theories_path, target_name=f"theory_{thid}")

    def download_eko(self, thid):
        """Download the EKO for a given theory ID"""
        thid = str(thid)
        remote = self.remote_ekos
        if thid not in remote:
            raise EkoNotFound(f"EKO for TheoryID {thid} is not available in the remote server")
        # Check that we have the theory we need
        target_path = self._ekos_path / f"eko_{int(thid)}.tar"
        download_file(remote[thid], target_path)

    def download_vp_output_file(self, filename, **kwargs):
        try:
            root_url = self.nnprofile['reports_root_url']
        except KeyError as e:
            raise LoadFailedError('Key report_root_url not found in nnprofile')
        try:
            url = root_url + filename
        except Exception as e:
            raise LoadFailedError(e) from e
        try:
            filename = pathlib.Path(filename)

            download_file(url, self._vp_cache() / filename, make_parents=True)
        except requests.HTTPError as e:
            if e.response.status_code == requests.codes.not_found:
                raise RemoteLoaderError(
                    f"Resource {filename} could not " f"be found on the validphys " f"server {url}"
                ) from e
            elif e.response.status_code == requests.codes.unauthorized:
                log.error(
                    "Could not access the validphys reports page "
                    "because the authentification is not provided. "
                    "Please, update your ~/.netrc file to contain the "
                    "following:\n\n"
                    f"machine {urls.urlsplit(root_url).netloc}\n"
                    f"    login nnpdf\n"
                    f"    password <PASSWORD>\n"
                )
            raise


class FallbackLoader(Loader, RemoteLoader):
    """A loader that first tries to find resources locally
    (calling Loader.check_*) and if it fails, it tries to download them
    (calling RemoteLoader.download_*)."""

    def make_checker(self, resource):
        # We are intercepting the check_
        orig = super().__getattribute__('check_' + resource)
        download = getattr(self, 'download_' + resource)

        @functools.wraps(orig)
        def f(*args, **kwargs):
            try:
                return orig(*args, **kwargs)
            except LoadFailedError as e:
                saved_exception = e
                log.info(
                    "Could not find a resource "
                    f"({resource}): {saved_exception}. "
                    f"Attempting to download it."
                )
                try:
                    download(*args, **kwargs)
                except RemoteLoaderError as e:
                    log.error(f"Failed to download resource: {e}")
                    raise e
                except LoadFailedError as e:
                    log.error(f"Resource not in the remote repository: {e}")
                    raise saved_exception
                except requests.RequestException as e:
                    log.error(f"There was a problem with the connection: {e}")
                    raise saved_exception from e

                except Exception as e:
                    # Simply raise these for now so we can find and fix them
                    raise e
                else:
                    return orig(*args, **kwargs)

        return f

    def __getattribute__(self, attr):
        token = 'check_'
        if attr.startswith(token):
            resname = attr[len(token) :]
            if hasattr(RemoteLoader, 'download_' + resname):
                return super().__getattribute__('make_checker')(resname)
        return super().__getattribute__(attr)
