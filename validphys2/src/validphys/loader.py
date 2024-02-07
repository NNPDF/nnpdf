# -*- coding: utf-8 -*-
"""
Resolve paths to useful objects, and query the existence of different resources
within the specified paths.
"""
import functools
from functools import cached_property
import logging
import mimetypes
import os
import os.path as osp
import pathlib
import pkgutil
import re
import shutil
import sys
import tempfile
from typing import List
import urllib.parse as urls

import requests

from reportengine import filefinder
from reportengine.compat import yaml
from validphys import lhaindex, pineparser
from validphys.commondataparser import parse_new_metadata
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
    peek_commondata_metadata,
)
from validphys.datafiles import path_vpdata
from validphys.utils import tempfile_cleaner

log = logging.getLogger(__name__)
NNPDF_DIR = "NNPDF"


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


def _get_nnpdf_profile(profile_path=None):
    """Returns the NNPDF profile as a dictionary

    If no ``profile_path`` is provided it will be autodiscovered in the following order:

    1. Environment variable $NNPDF_PROFILE_PATH
    2. ${XDG_CONFIG_HOME}/NNPDF/nnprofile.yaml (usually ~/.config/nnprofile)

    Any value not filled by 1 or 2 will then be filled by the default values
    found within the validphys python package `nnporfile_default.yaml`

    If ``nnpdf_share`` is set to the special key ``RELATIVE_TO_PYTHON``
    the python prefix (``Path(sys.prefix)/"share"/"NNPDF"``) will be used

    """
    yaml_reader = yaml.YAML(typ='safe', pure=True)

    home_config = pathlib.Path().home() / ".config"
    config_folder = pathlib.Path(os.environ.get("XDG_CONFIG_HOME", home_config)) / NNPDF_DIR

    # Set all default values
    profile_content = pkgutil.get_data("validphys", "nnprofile_default.yaml")
    profile_dict = yaml_reader.load(profile_content)
    # including the data_path to the validphys package
    profile_dict.setdefault("data_path", path_vpdata)

    # Look at profile path
    if profile_path is None:
        profile_path = os.environ.get("NNPDF_PROFILE_PATH", profile_path)

    # If profile_path is still none and there is a .config/NNPDF/nnprofile.yaml, read that
    if profile_path is None:
        if (config_nnprofile := config_folder / "nnprofile.yaml").exists():
            profile_path = config_nnprofile
        elif (config_nnprofile := config_folder / "nnprofile.yml").exists():
            profile_path = config_nnprofile

    if profile_path is not None:
        with open(profile_path, "r", encoding="utf-8") as f:
            profile_entries = yaml_reader.load(f)
            if profile_entries is not None:
                profile_dict.update(profile_entries)

    nnpdf_share = profile_dict.get("nnpdf_share")
    if nnpdf_share is None:
        if profile_path is not None:
            raise ValueError(
                f"`nnpdf_share` is not set in {profile_path}, please set it, e.g.: nnpdf_share: `.local/share/NNPDF`"
            )
        raise ValueError(
            "`nnpdf_share` not found in validphys, something is very wrong with the installation"
        )

    if nnpdf_share == "RELATIVE_TO_PYTHON":
        nnpdf_share = pathlib.Path(sys.prefix) / "share" / NNPDF_DIR

    # At this point nnpdf_share needs to be a path to somewhere
    nnpdf_share = pathlib.Path(nnpdf_share)

    # Make sure that we expand any ~ or ~<username>
    nnpdf_share = nnpdf_share.expanduser()

    # Make sure we can either write to this directory or it exists
    try:
        nnpdf_share.mkdir(exist_ok=True, parents=True)
    except PermissionError as e:
        raise FileNotFoundError(
            f"{nnpdf_share} does not exist and you haven't got permissions to create it!"
        ) from e

    # Now read all paths and define them as relative to nnpdf_share (unless given as absolute)
    for var in ["results_path", "theories_path", "validphys_cache_path", "hyperscan_path"]:
        # if there are any problems setting or getting these variable erroring out is more than justified
        absolute_var = nnpdf_share / pathlib.Path(profile_dict[var]).expanduser()
        profile_dict[var] = absolute_var.absolute().as_posix()

    return profile_dict


class LoaderBase:
    """
    Base class for the NNPDF loader.
    It can take as input a profile dictionary from which all data can be read.
    It is possible to override the datapath and resultpath when the class is instantiated.
    """

    def __init__(self, profile=None):
        if not isinstance(profile, dict):
            # If profile is a path, a str or None, read it from the default path
            profile = _get_nnpdf_profile(profile)

        # Retrieve important paths from the profile if not given
        datapath = pathlib.Path(profile["data_path"])
        theories_path = pathlib.Path(profile["theories_path"])
        resultspath = pathlib.Path(profile["results_path"])

        if not datapath.exists():
            raise LoaderError(f"The data path {datapath} does not exist.")

        # Create the theories and results paths if they don't exist already
        theories_path.mkdir(exist_ok=True, parents=True)
        resultspath.mkdir(exist_ok=True, parents=True)

        # And save them up
        self.datapath = datapath
        self._theories_path = theories_path
        self.resultspath = resultspath
        self._old_commondata_fits = set()
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


def rebuild_commondata_without_cuts(filename_with_cuts, cuts, datapath_filename, newpath):
    """Take a CommonData file that is stored with the cuts applied
    and write another file with no cuts. The points that were not present in
    the original file have the same kinematics as the file in
    ``datapath_filename``, which must correspond to the original CommonData
    file which does not have the cuts applied. However, to avoid confusion, the
    values and uncertainties are all set to zero. The new file is written
    to ``newpath``.
    """

    metadata = peek_commondata_metadata(datapath_filename)
    if cuts is None:
        shutil.copy2(filename_with_cuts, newpath)
        return

    index_pattern = re.compile(r'(?P<startspace>\s*)(?P<index>\d+)')
    data_line_pattern = re.compile(
        r'\s*(?P<index>\d+)' r'\s+(?P<process_type>\S+)\s+' r'(?P<kinematics>(\s*\S+){3})\s+'
    )
    mask = cuts.load()
    maskiter = iter(mask)
    ndata = metadata.ndata
    nsys = metadata.nsys

    next_index = next(maskiter)
    with open(filename_with_cuts, 'r') as fitfile, open(datapath_filename) as dtfile, open(
        newpath, 'w'
    ) as newfile:
        newfile.write(dtfile.readline())
        # discard this line
        fitfile.readline()
        for i in range(1, ndata + 1):
            # You gotta love mismatched indexing
            if i - 1 == next_index:
                line = fitfile.readline()
                line = re.sub(index_pattern, rf'\g<startspace>{i}', line, count=1)
                newfile.write(line)
                next_index = next(maskiter, None)
                # drop the data file line
                dtfile.readline()
            else:
                line = dtfile.readline()
                # check that we know where we are
                m = re.match(index_pattern, line)
                assert int(m.group('index')) == i
                # We have index, process type, and 3*kinematics
                # that we would like to keep.
                m = re.match(data_line_pattern, line)
                newfile.write(line[: m.end()])
                # And value, stat, *sys that we want to drop
                # Do not use string join to keep up with the ugly format
                # This should really be nan's, but the c++ streams that could read this
                # do not have the right interface.
                # https://stackoverflow.com/questions/11420263/is-it-possible-to-read-infinity-or-nan-values-using-input-streams
                zeros = '-0\t' * (2 + 2 * nsys)
                newfile.write(f'{zeros}\n')


# TODO: Deprecate get methods?
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
    @functools.lru_cache()
    def available_theories(self):
        """Return a string token for each of the available theories"""
        theory_token = 'theory_'
        return {
            folder.name[len(theory_token) :]
            for folder in self._theories_path.glob(theory_token + '*')
        }

    @property
    @functools.lru_cache()
    def available_datasets(self):
        data_str = "DATA_"
        # We filter out the positivity sets here
        return {
            file.stem[len(data_str) :]
            for file in self.commondata_folder.glob(f'{data_str}*.dat')
            if not file.stem.startswith((f"{data_str}POS", f"{data_str}INTEG"))
        }

    @property
    @functools.lru_cache()
    def available_pdfs(self):
        return lhaindex.expand_local_names('*')

    @property
    def commondata_folder(self):
        return self.datapath / 'commondata'

    def check_commondata(
        self, setname, sysnum=None, use_fitcommondata=False, fit=None, variant=None
    ):
        if use_fitcommondata:
            if not fit:
                raise LoadFailedError("Must specify a fit when setting use_fitcommondata")
            datafilefolder = (fit.path / 'filter') / setname
            newpath = datafilefolder / f'FILTER_{setname}.dat'
            if not newpath.exists():
                oldpath = datafilefolder / f'DATA_{setname}.dat'
                if not oldpath.exists():
                    raise DataNotFoundError(
                        f"Either {newpath} or {oldpath} are needed with `use_fitcommondata`"
                    )
                # This is to not repeat all the error handling stuff
                basedata = self.check_commondata(setname, sysnum=sysnum)
                basedata_path = basedata.datafile
                cuts = self.check_fit_cuts(basedata, fit=fit)

                if fit not in self._old_commondata_fits:
                    self._old_commondata_fits.add(fit)
                    log.warning(
                        f"Found fit using old commondata export settings: "
                        f"'{fit}'. The commondata that are used in this run "
                        "will be updated now."
                        "Please consider re-uploading it."
                    )
                    log.warning("Points that do not pass the cuts are set to zero!")

                log.info(f"Upgrading filtered commondata. Writing {newpath}")
                rebuild_commondata_without_cuts(oldpath, cuts, basedata_path, newpath)
            datafile = newpath
        else:
            datafile = self.commondata_folder / f'DATA_{setname}.dat'

        if not datafile.exists():
            # TODO: if not old data found, maybe thi sis a new data
            # The new data goes into folder inside `self.commondata_folder`
            # this usually corresponds to <validphys_code>/datafiles/commondata
            setfolder, observable_name = setname.rsplit("_", 1)

            # TODO
            if not self.commondata_folder.with_name("new_commondata").exists():
                raise DataNotFoundError("new_commondata folder missing in this branch!")

            metadata_file = (
                self.commondata_folder.with_name("new_commondata") / setfolder / "metadata.yaml"
            )

            # If the metadata file doesn't exist either, then error out
            if not metadata_file.exists():
                raise DataNotFoundError(
                    f"""The CommonData set {setname} could not be found
as old ({datafile})
or new ({metadata_file})"""
                )

            # Get the instance of ObservableMetaData
            metadata = parse_new_metadata(metadata_file, observable_name, variant=variant)

            return CommonDataSpec(None, None, None, name=setname, metadata=metadata, legacy=False)

        if sysnum is None:
            sysnum = 'DEFAULT'
        sysfile = self.commondata_folder / 'systypes' / ('SYSTYPE_%s_%s.dat' % (setname, sysnum))

        if not sysfile.exists():
            raise SysNotFoundError(
                "Could not find systype %s for dataset '%s'. File %s does not exist."
                % (sysnum, setname, sysfile)
            )

        plotfiles = []

        metadata = peek_commondata_metadata(datafile)
        process_plotting_root = self.commondata_folder / f'PLOTTINGTYPE_{metadata.process_type}'
        type_plotting = (
            process_plotting_root.with_suffix('.yml'),
            process_plotting_root.with_suffix('.yaml'),
        )

        data_plotting_root = self.commondata_folder / f'PLOTTING_{setname}'

        data_plotting = (
            data_plotting_root.with_suffix('.yml'),
            data_plotting_root.with_suffix('.yaml'),
        )
        # TODO: What do we do when both .yml and .yaml exist?
        for tp in (type_plotting, data_plotting):
            for p in tp:
                if p.exists():
                    plotfiles.append(p)
        if setname != metadata.name:
            raise InconsistentMetaDataError(
                f"The name found in the CommonData file, {metadata.name}, did "
                f"not match the dataset name, {setname}."
            )
        return CommonDataSpec(datafile, sysfile, plotfiles, name=setname, metadata=metadata)

    @functools.lru_cache()
    def check_theoryID(self, theoryID):
        theoryID = str(theoryID)
        theopath = self._theories_path / f"theory_{theoryID}"
        if not theopath.exists():
            raise TheoryNotFound(
                "Could not find theory %s. Folder '%s' not found" % (theoryID, theopath)
            )
        return TheoryIDSpec(theoryID, theopath, self.theorydb_file)

    @property
    def theorydb_file(self):
        """Checks theory db file exists and returns path to it"""
        dbpath = self.datapath / 'theory.db'
        if not dbpath.is_file():
            raise TheoryDataBaseNotFound(f"could not find theory.db. File not found at {dbpath}")
        return dbpath

    def get_commondata(self, setname, sysnum):
        """Get a Commondata from the set name and number."""
        # TODO: check where this is used
        # as this might ignore cfactors or variants
        cd = self.check_commondata(setname, sysnum)
        return cd.load()

    #   @functools.lru_cache()
    def check_fktable(self, theoryID, setname, cfac):
        _, theopath = self.check_theoryID(theoryID)
        fkpath = theopath / 'fastkernel' / ('FK_%s.dat' % setname)
        if not fkpath.exists():
            raise FKTableNotFound(
                "Could not find FKTable for set '%s'. File '%s' not found" % (setname, fkpath)
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

    def check_fkyaml(self, name, theoryID, cfac):
        """Load a pineappl fktable in the old commondata format
        Receives a yaml file describing the fktables necessary for a given observable
        the theory ID and the corresponding cfactors.
        The cfactors should correspond directly to the fktables, the "compound folder"
        is not supported for pineappl theories. As such, the name of the cfactor is expected to be
            CF_{cfactor_name}_{fktable_name}
        """
        theory = self.check_theoryID(theoryID)
        if (theory.path / "compound").exists():
            raise LoadFailedError(f"New theories (id=${theoryID}) do not accept compound files")

        fkpath = (theory.yamldb_path / name).with_suffix(".yaml")
        theory_metadata, _ = pineparser.get_yaml_information(fkpath, theory.path)
        return self.check_fk_from_theory_metadata(theory_metadata, theoryID, cfac)

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
        data = yaml.safe_load(yaml_format)
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

    def get_fktable(self, theoryID, setname, cfac):
        fkspec = self.check_fktable(theoryID, setname, cfac)
        return fkspec.load()

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
        if cd.legacy:
            if th.is_pineappl():
                fk, _ = self.check_fkyaml(setname, theoryID, [])
            else:
                fk = self.check_fktable(theoryID, setname, [])
        else:
            fk, _ = self.check_fk_from_theory_metadata(cd.metadata.theory, theoryID)
        return cd, fk, th

    def check_posset(self, theoryID, setname, postlambda):
        """Load a positivity dataset"""
        cd, fk, th = self._check_lagrange_multiplier_set(theoryID, setname)
        return PositivitySetSpec(setname, cd, fk, postlambda, th)

    def check_integset(self, theoryID, setname, postlambda):
        """Load an integrability dataset"""
        cd, fk, th = self._check_lagrange_multiplier_set(theoryID, setname)
        return IntegrabilitySetSpec(setname, cd, fk, postlambda, th)

    def get_posset(self, theoryID, setname, postlambda):
        return self.check_posset(theoryID, setname, postlambda).load()

    def check_fit(self, fitname):
        resultspath = self.resultspath
        if fitname != osp.basename(fitname):
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
        if hyperscan_name != osp.basename(hyperscan_name):
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
        return [
            Rule(inp, defaults=defaults, theory_parameters=th_params, loader=self)
            for inp in default_filter_rules_input()
        ]

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

        theoryno, _ = theoryid

        # TODO:
        # The dataset is checked twice, once here
        # and once by config in produce_commondata
        # once of the two __must__ be superfluous
        # note that both use information from dataset_input
        commondata = self.check_commondata(
            name, sysnum, use_fitcommondata=use_fitcommondata, fit=fit, variant=variant
        )

        if commondata.legacy:
            if theoryid.is_pineappl():
                # If it is a pineappl theory, use the pineappl reader
                fkspec, op = self.check_fkyaml(name, theoryno, cfac)
            else:
                try:
                    fkspec, op = self.check_compound(theoryno, name, cfac)
                except CompoundNotFound:
                    fkspec = self.check_fktable(theoryno, name, cfac)
                    op = None
        else:
            # New commondata files work _only_ with pineappl theory
            if not theoryid.is_pineappl():
                raise ValueError(
                    f"New commondata files accept only pineappl theories (used:{theoryid.id})"
                )

            if (thmeta := commondata.metadata.theory) is None:
                raise TheoryMetadataNotFound(f"No theory metadata found for {name}")

            fkspec, op = self.check_fk_from_theory_metadata(thmeta, theoryno, cfac)

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
        )

    def check_experiment(self, name: str, datasets: List[DataSetSpec]) -> DataGroupSpec:
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
        >>> ds = l.check_dataset("NMC", theoryid=53, cuts="internal")
        >>> exp = l.check_experiment("My DataGroupSpec Name", [ds])
        """
        if not isinstance(datasets, list):
            raise TypeError("Must specify a list of DataSetSpec objects to use")

        return DataGroupSpec(name, datasets)

    def check_pdf(self, name):
        if lhaindex.isinstalled(name):
            return PDF(name)
        raise PDFNotFound(name)

    def get_pdf(self, name):
        return self.check_pdf(name).load()

    def check_fit_cuts(self, commondata, fit):
        setname = commondata.name
        if fit is None:
            raise TypeError("Must specify a fit to use the cuts.")
        if not isinstance(fit, FitSpec):
            fit = self.check_fit(fit)
        _, fitpath = fit
        p = (fitpath / 'filter') / setname / ('FKMASK_' + setname + '.dat')
        if not p.parent.exists():
            raise CutsNotFound(f"Bad filter configuration. Could not find {p.parent}")
        if not p.exists():
            p = None
        return Cuts(commondata, p)

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


def download_file(url, stream_or_path, make_parents=False):
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

        download_target = tempfile.NamedTemporaryFile(
            delete=False, dir=p.parent, prefix=p.name, suffix='.part'
        )

        with download_target as f:
            _download_and_show(response, f)
        shutil.move(download_target.name, p)
    else:
        log.info("Downloading %s.", url)
        _download_and_show(response, stream_or_path)


def download_and_extract(url, local_path):
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
        shutil.unpack_archive(t.name, extract_dir=local_path)
    except:
        log.error(
            f"The original archive at {t.name} was only extracted partially at \n{local_path}"
        )
        raise
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
            raise RemoteLoaderError(
                "Failed to fetch remote %s index %s: %s" % (thing, index_url, e)
            ) from e

        try:
            info = resp.json()['files']
        except Exception as e:
            raise RemoteLoaderError(
                "Malformed index %s. Expecting json with a key 'files': %s" % (index_url, e)
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
    @functools.lru_cache()
    def remote_fits(self):
        return self.remote_files(self.fit_urls, self.fit_index, thing="fits")

    @property
    @functools.lru_cache()
    def remote_hyperscans(self):
        return self.remote_files(self.hyperscan_url, self.hyperscan_index, thing="hyperscan")

    @property
    @functools.lru_cache()
    def remote_theories(self):
        token = 'theory_'
        rt = self.remote_files(self.theory_urls, self.theory_index, thing="theories")
        return {k[len(token) :]: v for k, v in rt.items()}

    @property
    @functools.lru_cache()
    def remote_nnpdf_pdfs(self):
        return self.remote_files(self.nnpdf_pdfs_urls, self.nnpdf_pdfs_index, thing="PDFs")

    @cached_property
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
    def lhapdf_pdfs(self):
        return lhaindex.expand_index_names('*')

    @property
    def nnpdf_pdfs(self):
        return list(self.remote_nnpdf_pdfs)

    @property
    def downloadable_pdfs(self):
        return set((*self.lhapdf_pdfs, *self.downloadable_fits, *self.nnpdf_pdfs))

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
        thid = str(thid)
        remote = self.remote_theories
        if thid not in remote:
            raise TheoryNotFound("Theory %s not available." % thid)
        download_and_extract(remote[thid], self._theories_path)

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
