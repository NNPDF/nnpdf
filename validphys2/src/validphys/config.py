# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:43:10 2016

@author: Zahari Kassabov
"""
import logging
import pathlib
import functools
import inspect
import numbers
import copy
import os
from importlib.resources import read_text, contents

from collections import ChainMap, defaultdict
from collections.abc import Mapping, Sequence

from reportengine import configparser
from reportengine.environment import Environment, EnvironmentError_
from reportengine.configparser import (
    ConfigError,
    element_of,
    _parse_func,
    record_from_defaults,
)
from reportengine.helputils import get_parser_type
from reportengine.namespaces import NSList
from reportengine import report
from reportengine.compat import yaml

from validphys.core import (
    DataGroupSpec,
    DataSetInput,
    ExperimentInput,
    CutsPolicy,
    MatchedCuts,
    SimilarCuts,
    ThCovMatSpec,
)
from validphys.loader import (
    Loader,
    LoaderError,
    LoadFailedError,
    DataNotFoundError,
    PDFNotFound,
    FallbackLoader,
    InconsistentMetaDataError,
)
from validphys.gridvalues import LUMI_CHANNELS

from validphys.paramfits.config import ParamfitsConfig

from validphys.theorycovariance.theorycovarianceutils import process_lookup
from validphys.plotoptions import get_info

import validphys.scalevariations

log = logging.getLogger(__name__)


class Environment(Environment):
    """Container for information to be filled at run time"""

    def __init__(
        self,
        *,
        datapath=None,
        resultspath=None,
        this_folder=None,
        net=True,
        upload=False,
        **kwargs,
    ):
        if this_folder:
            self.this_folder = pathlib.Path(this_folder)

        if net:
            loader_class = FallbackLoader
        else:
            loader_class = Loader

        try:
            self.loader = loader_class(datapath, resultspath)
        except LoaderError as e:
            log.error(
                "Failed to find the paths. These are configured "
                "in the nnprofile settings or with the "
                "--datapath and --resultpath options."
            )
            raise EnvironmentError_(e) from e
        self.deta_path = self.loader.datapath
        self.results_path = self.loader.resultspath

        self.upload = upload
        super().__init__(**kwargs)


def _id_with_label(f):
    f = _parse_func(f)

    def parse_func(self, item, **kwargs):
        if not isinstance(item, dict):
            return f(self, item, **kwargs)
        keydiff = item.keys() - {"id", "label"}

        if keydiff or not "id" in item:
            unrecognized = f" Unrecognized: {keydiff}" if keydiff else ""
            raise ConfigError(
                f"'{item}' must be a single id, or a mapping "
                f"with keys 'id', 'label.{unrecognized}'"
            )
        id = item["id"]
        val = f(self, id, **kwargs)
        if "label" in item:
            val.label = str(item["label"])
        return val

    currsig = inspect.signature(parse_func)
    origsig = inspect.signature(f)
    parse_func = functools.wraps(f)(parse_func)

    params = [
        *list(currsig.parameters.values())[:2],
        *list(origsig.parameters.values())[2:],
    ]

    parse_func.__signature__ = inspect.Signature(parameters=params)

    labeldoc = (" Either just an id %s, or a mapping " "with 'id' and 'label'.") % (
        get_parser_type(f),
    )
    if parse_func.__doc__ is None:
        parse_func.__doc__ = labeldoc
    else:
        parse_func.__doc__ += labeldoc

    return parse_func


class CoreConfig(configparser.Config):
    @property
    def loader(self):
        return self.environment.loader

    @element_of("pdfs")
    @_id_with_label
    def parse_pdf(self, name: str):
        """A PDF set installed in LHAPDF."""
        try:
            pdf = self.loader.check_pdf(name)
        except PDFNotFound as e:
            raise ConfigError(
                "Bad PDF: {} not installed".format(name),
                name,
                self.loader.available_pdfs,
            ) from e
        except LoaderError as e:
            raise ConfigError(e) from e

        # Check that we know how to compute errors
        try:
            pdf.nnpdf_error
        except NotImplementedError as e:
            raise ConfigError(str(e))
        return pdf

    @element_of("theoryids")
    @_id_with_label
    def parse_theoryid(self, theoryID: (str, int)):
        """A number corresponding to the database theory ID where the
        corresponding theory folder is installed in te data directory."""
        try:
            return self.loader.check_theoryID(theoryID)
        except LoaderError as e:
            raise ConfigError(
                str(e),
                theoryID,
                self.loader.available_theories,
                display_alternatives="all",
            )

    def parse_use_cuts(self, use_cuts: (bool, str)):
        """Whether to filter the points based on the cuts applied in the fit,
        or the whole data in the dataset. The possible options are:

        - internal: Calculate the cuts based on the existing rules. This is
          the default.

        - fromfit: Read the cuts stored in the fit.

        - nocuts: Use the whole dataset.
        """
        # The lower is an aesthetic preference...
        valid_cuts = {c.value for c in CutsPolicy}
        if isinstance(use_cuts, bool):
            if use_cuts:
                res = CutsPolicy.FROMFIT
            else:
                res = CutsPolicy.NOCUTS
            log.warning(
                "Setting a boolean for `use_cuts` is deprecated. "
                f"The available values are {valid_cuts} and the default "
                f"value is 'internal'. Your input ('{use_cuts}') is "
                f"equivalent to '{res}'."
            )
        elif isinstance(use_cuts, str) and use_cuts in valid_cuts:
            res = CutsPolicy(use_cuts)
        else:
            raise ConfigError(
                f"Invalid use_cuts setting: '{use_cuts}'.", use_cuts, valid_cuts
            )

        return res

    # TODO: load fit config from here
    @element_of("fits")
    @_id_with_label
    def parse_fit(self, fit: str):
        """A fit in the results folder, containing at least a valid filter result."""
        try:
            return self.loader.check_fit(fit)
        except LoadFailedError as e:
            raise ConfigError(str(e), fit, self.loader.available_fits)

    def produce_fitcontext(self, fitinputcontext, fitpdf):
        """Set PDF, theory ID and data input from the fit config"""

        return dict(**fitinputcontext, **fitpdf)

    def produce_fitinputcontext(self, fit):
        """Like ``fitcontext`` but without setting the PDF"""
        _, theory = self.parse_from_("fit", "theory", write=False)
        thid = theory["theoryid"]

        data_input = self._parse_data_input_from_(
            "fit", {"theoryid": thid}
        )
        return {"theoryid": thid, "data_input": data_input}

    def produce_fitpdf(self, fit):
        """Like ``fitcontext`` only setting the PDF"""
        with self.set_context(ns=self._curr_ns.new_child({"fit": fit})):
            _, pdf = self.parse_from_("fit", "pdf", write=False)
        return {"pdf": pdf}

    def produce_fitunderlyinglaw(self, fit):
        """Reads closuretest: fakepdf from fit config file and passes as
        pdf
        """
        with self.set_context(ns=self._curr_ns.new_child({"fit": fit})):
            _, datacuts = self.parse_from_("fit", "closuretest", write=False)
        underlyinglaw = self.parse_pdf(datacuts["fakepdf"])
        return {"pdf": underlyinglaw}

    def produce_multiclosure_underlyinglaw(self, fits):
        """Produce the underlying law for a set of fits. This allows a single t0
        like covariance matrix to be loaded for all fits, for use with
        statistical estimators on multiple closure fits. If the fits don't all
        have the same underlying law then an error is raised, offending fit is
        identified.
        """
        # could use comprehension here but more useful to find offending fit
        laws = set()
        for fit in fits:
            try:
                closuretest_spec = fit.as_input()["closuretest"]
            except KeyError as e:
                raise ConfigError(
                    f"fit: {fit} does not have a `closuretest` namespace in " "runcard"
                ) from e
            try:
                laws.add(closuretest_spec["fakepdf"])
            except KeyError as e:
                raise ConfigError(
                    f"fit: {fit} does not have `fakepdf` specified in the "
                    "closuretest namespace in runcard."
                ) from e

        if len(laws) != 1:
            raise ConfigError(
                "Did not find unique underlying law from fits, "
                f"instead found: {laws}"
            )
        return self.parse_pdf(laws.pop())


    def produce_basisfromfit(self, fit):
        """Set the basis from fit config. In the fit config file the basis
        is set using the key ``fitbasis``, but it is exposed to validphys
        as ``basis``.

        The name of this production rule is intentionally
        set to not conflict with the existing ``fitbasis`` runcard key.

        """
        with self.set_context(ns=self._curr_ns.new_child({"fit": fit})):
            _, fitting = self.parse_from_("fit", "fitting", write=False)
        basis = fitting["fitbasis"]
        return {"basis": basis}


    def produce_fitpdfandbasis(self, fitpdf, basisfromfit):
        """ Set the PDF and basis from the fit config. """
        return {**fitpdf, **basisfromfit}


    @element_of("dataset_inputs")
    def parse_dataset_input(self, dataset: Mapping):
        """The mapping that corresponds to the dataset specifications in the
        fit files"""
        known_keys = {"dataset", "sys", "cfac", "frac", "weight"}
        try:
            name = dataset["dataset"]
            if not isinstance(name, str):
                raise ConfigError(f"'dataset' must be a string, not {type(name)}")
        except KeyError:
            raise ConfigError(
                "'dataset' must be a mapping with " "'dataset' and 'sysnum'"
            )

        sysnum = dataset.get("sys")
        cfac = dataset.get("cfac", tuple())
        frac = dataset.get("frac", 1)
        if not isinstance(frac, numbers.Real):
            raise ConfigError(f"'frac' must be a number, not '{frac}'")
        if frac < 0 or frac > 1:
            raise ConfigError(f"'frac' must be between 0 and 1 not '{frac}'")
        weight = dataset.get("weight", 1)
        if not isinstance(weight, numbers.Real):
            raise ConfigError(f"'weight' must be a number, not '{weight}'")
        if weight < 0:
            raise ConfigError(f"'weight' must be greater than zero not '{weight}'")
        kdiff = dataset.keys() - known_keys
        for k in kdiff:
            # Abuse ConfigError to get the suggestions.
            log.warning(
                ConfigError(f"Key '{k}' in dataset_input not known.", k, known_keys)
            )
        return DataSetInput(name=name, sys=sysnum, cfac=cfac, frac=frac, weight=weight)

    def parse_use_fitcommondata(self, do_use: bool):
        """Use the commondata files in the fit instead of those in the data
        directory."""
        return do_use

    def produce_commondata(self, *, dataset_input, use_fitcommondata=False, fit=None):
        """Produce a CommondataSpec from a dataset input"""

        name = dataset_input.name
        sysnum = dataset_input.sys
        try:
            return self.loader.check_commondata(
                setname=name,
                sysnum=sysnum,
                use_fitcommondata=use_fitcommondata,
                fit=fit,
            )
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets) from e
        except LoadFailedError as e:
            raise ConfigError(e) from e
        except InconsistentMetaDataError as e:
            raise ConfigError(e) from e

    def parse_cut_similarity_threshold(self, th: numbers.Real):
        """Maximum relative ratio when using `fromsimilarpredictons` cuts."""
        return th

    def produce_cuts(
        self, *, commondata, use_cuts, rules, fit=None, theoryid=None,
    ):
        """Obtain cuts for a given dataset input, based on the
        appropriate policy."""
        # TODO: Put this bit of logic into loader.check_cuts
        if use_cuts is CutsPolicy.NOCUTS:
            return None
        elif use_cuts is CutsPolicy.FROMFIT:
            if not fit:
                raise ConfigError(
                    "Setting 'use_cuts' to 'fromfit' requires "
                    "specifying a fit on which filter "
                    "has been executed, e.g.\nfit : NNPDF30_nlo_as_0118"
                )
            name = commondata.name
            try:
                return self.loader.check_fit_cuts(name, fit)
            except LoadFailedError as e:
                raise ConfigError(e) from e
        elif use_cuts is CutsPolicy.INTERNAL:
            if not theoryid:
                raise ConfigError("theoryid must be specified for internal cuts")
            return self.loader.check_internal_cuts(commondata, rules)
        elif (
            use_cuts is CutsPolicy.FROM_CUT_INTERSECTION_NAMESPACE
            or use_cuts is CutsPolicy.FROM_SIMILAR_PREDICTIONS_NAMESPACE
        ):
            cut_list = []
            _, nss = self.parse_from_(None, "cuts_intersection_spec", write=False)
            self._check_dataspecs_type(nss)
            if not nss:
                raise ConfigError(
                    "'cuts_intersection_spec' must contain at least one namespace."
                )

            ns_cut_inputs = {"commondata": commondata, "use_cuts": CutsPolicy.INTERNAL}
            for ns in nss:
                with self.set_context(
                    ns=self._curr_ns.new_child({**ns, **ns_cut_inputs})
                ):
                    _, nscuts = self.parse_from_(None, "cuts", write=False)
                    cut_list.append(nscuts)
            ndata = commondata.ndata
            matched_cuts = MatchedCuts(cut_list, ndata=ndata)
            if use_cuts is CutsPolicy.FROM_CUT_INTERSECTION_NAMESPACE:
                return matched_cuts
            else:
                if len(nss) != 2:
                    raise ConfigError("Can only work with two namespaces")
                _, cut_similarity_threshold = self.parse_from_(
                    None, "cut_similarity_threshold", write=False
                )
                name = commondata.name
                inps = []
                for i, ns in enumerate(nss):
                    with self.set_context(ns=self._curr_ns.new_child({**ns,})):
                        # TODO: find a way to not duplicate this and use a dict
                        # instead of a linear search
                        _, dins = self.parse_from_(None, "dataset_inputs", write=False)
                    try:
                        di = next(d for d in dins if d.name == name)
                    except StopIteration as e:
                        raise ConfigError(
                            f"cuts_intersection_spec namespace {i}: dataset inputs must define {name}"
                        ) from e

                    with self.set_context(
                        ns=self._curr_ns.new_child(
                            {
                                "dataset_input": di,
                                "use_cuts": CutsPolicy.FROM_CUT_INTERSECTION_NAMESPACE,
                                "cuts": matched_cuts,
                                **ns,
                            }
                        )
                    ):
                        _, ds = self.parse_from_(None, "dataset", write=False)
                        _, pdf = self.parse_from_(None, "pdf", write=False)
                    inps.append((ds, pdf))
                return SimilarCuts(tuple(inps), cut_similarity_threshold)

        raise TypeError("Wrong use_cuts")

    def produce_dataset(
        self,
        *,
        rules,
        dataset_input,
        theoryid,
        cuts,
        use_fitcommondata=False,
        fit=None,
        check_plotting: bool = False,
    ):
        """Dataset specification from the theory and CommonData.
           Use the cuts from the fit, if provided. If check_plotting is set to
           True, attempt to lod and check the PLOTTING files
           (note this may cause a noticeable slowdown in general)."""
        name = dataset_input.name
        sysnum = dataset_input.sys
        cfac = dataset_input.cfac
        frac = dataset_input.frac
        weight = dataset_input.weight

        try:
            ds = self.loader.check_dataset(
                rules=rules,
                name=name,
                sysnum=sysnum,
                theoryid=theoryid,
                cfac=cfac,
                cuts=cuts,
                frac=frac,
                use_fitcommondata=use_fitcommondata,
                fit=fit,
                weight=weight,
            )
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets)

        except LoadFailedError as e:
            raise ConfigError(e)

        if check_plotting:
            from validphys.plotoptions import get_info

            # normalize=True should check for more stuff
            get_info(ds, normalize=True)
            if not ds.commondata.plotfiles:
                log.warning(f"Plotting files not found for: {ds}")
        return ds

    @configparser.element_of("experiments")
    def parse_experiment(self, experiment: dict):
        """A set of datasets where correlated systematics are taken
           into account. It is a mapping where the keys are the experiment
           name 'experiment' and a list of datasets."""
        try:
            name, datasets = experiment["experiment"], experiment["datasets"]
        except KeyError as e:
            raise ConfigError(
                "'experiment' must be a mapping with "
                "'experiment' and 'datasets', but %s is missing" % e
            ) from e

        dsinputs = [self.parse_dataset_input(ds) for ds in datasets]


        return self.produce_data(group_name=name, data_input=dsinputs)

    @configparser.element_of("experiment_inputs")
    def parse_experiment_input(self, ei: dict):
        """The mapping that corresponds to the experiment specification in the
        fit config files. Currently, this needs to be combined with
        ``experiment_from_input`` to yield an useful result."""
        try:
            name = ei["experiment"]
        except KeyError as e:
            raise ConfigError(f"experiment_input must have an 'experiment' key") from e

        try:
            datasets = ei["datasets"]
        except KeyError as e:
            raise ConfigError(f"experiment_input must have an 'datasets' key") from e

        return ExperimentInput(name=name, datasets=datasets)

    # TODO: Do away with the mapping and make the conversion implicitly
    def produce_experiment_from_input(
        self, experiment_input, theoryid, use_cuts, fit=None
    ):
        """Return a mapping containing a single experiment from an experiment
        input. NOTE: This might be deprecated in the future."""
        return {
            "experiment": self.parse_experiment(
                experiment_input.as_dict(),
                theoryid=theoryid,
                use_cuts=use_cuts,
                fit=fit,
            )
        }

    @configparser.explicit_node
    def produce_covariance_matrix(self, use_pdferr: bool = False):
        """Modifies which action is used as covariance_matrix depending on
        the flag `use_pdferr`
        """
        from validphys import covmats

        if use_pdferr:
            return covmats.pdferr_plus_covmat
        else:
            return covmats.covmat

    @configparser.explicit_node
    def produce_dataset_inputs_covariance_matrix(self, use_pdferr: bool = False):
        """Modifies which action is used as experiment_covariance_matrix
        depending on the flag `use_pdferr`
        """
        from validphys import covmats

        if use_pdferr:
            return covmats.pdferr_plus_dataset_inputs_covmat
        else:
            return covmats.dataset_inputs_covmat

    # TODO: Do this better and elsewhere
    @staticmethod
    def _check_dataspecs_type(dataspecs):
        if not isinstance(dataspecs, Sequence):
            raise ConfigError(
                "dataspecs should be a sequence of mappings, not "
                f"{type(dataspecs).__name__}"
            )

        for spec in dataspecs:
            if not isinstance(spec, Mapping):
                raise ConfigError(
                    "dataspecs should be a sequence of mappings, "
                    f" but {spec} is {type(spec).__name__}"
                )

    def produce_matched_datasets_from_dataspecs(self, dataspecs):
        """Take an arbitrary list of mappings called dataspecs and
        return a new list of mappings called dataspecs constructed as follows.

        From each of the original dataspecs, resolve the key `process`, and
        all the experiments and datasets therein.

        Compute the intersection of the dataset names, and for each element in
        the intersection construct a mapping with the follwing keys:

            - process : A string with the common process name.
            - experiment_name : A string with the common experiment name.
            - dataset_name : A string with the common dataset name.
            - dataspecs : A list of mappinngs matching the original
              "dataspecs". Each mapping contains:
                * dataset: A dataset with the name data_set name and the
                properties (cuts, theory, etc) corresponding to the original
                dataspec.
                * dataset_input: The input line used to build dataset.
                * All the other keys in the original dataspec.
        """
        self._check_dataspecs_type(dataspecs)
        all_names = []
        for spec in dataspecs:

            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, data_input = self.parse_from_(None, "data_input", write=False)
                names = {
                    (process_lookup(dsin.name), dsin.name): dsin for dsin in data_input
                }

                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))
        res = []
        for k in used_set:
            inres = {"process": k[0], "dataset_name": k[1]}
            # TODO: Should this have the same name?
            inner_spec_list = inres["dataspecs"] = []
            for ispec, spec in enumerate(dataspecs):
                # Passing spec by referene
                d = ChainMap({"dataset_input": all_names[ispec][k],}, spec)
                inner_spec_list.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x["process"], x["dataset_name"]))
        return res

    def produce_matched_positivity_from_dataspecs(self, dataspecs):
        """Like produce_matched_datasets_from_dataspecs but for positivity datasets."""
        self._check_dataspecs_type(dataspecs)
        all_names = []
        for spec in dataspecs:
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, pos = self.parse_from_(None, "posdatasets", write=False)
                names = {(p.name): (p) for p in pos}
                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))

        res = []
        for k in used_set:
            inres = {"posdataset_name": k}
            # TODO: Should this have the same name?
            l = inres["dataspecs"] = []
            for ispec, spec in enumerate(dataspecs):
                # Passing spec by referene
                d = ChainMap({"posdataset": all_names[ispec][k],}, spec)
                l.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x["posdataset_name"]))
        return res

    def produce_dataspecs_with_matched_cuts(self, dataspecs):
        """Take a list of namespaces (dataspecs), resolve ``dataset`` within
        each of them, and return another list of dataspecs where the datasets
        all have the same cuts, corresponding to the intersection of the
        selected points. All the datasets must have the same name (i.e.
        correspond with the same experimental measurement), but can otherwise
        differ, for example in the theory used for the experimental
        predictions.

        This rule can be combined with ``matched_datasets_from_dataspecs``.
        """
        self._check_dataspecs_type(dataspecs)
        if not dataspecs:
            return dataspecs
        # Can now assume we have at least one element
        cutlist = []
        dslist = []
        names = set()
        for spec in dataspecs:
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, ds = self.parse_from_(None, "dataset", write=False)
            dslist.append(ds)
            cutlist.append(ds.cuts)
            names.add(ds.name)

        lnames = len(names)
        if lnames != 1:
            raise ConfigError(
                "Each dataspec must have a dataset with the same"
                f"name, but got {lnames} different ones: {names}"
            )

        ndata = ds.commondata.ndata
        matched_cuts = MatchedCuts(cutlist, ndata=ndata)
        res = []
        for spec, ds in zip(dataspecs, dslist):
            newds = copy.copy(ds)
            newds.cuts = matched_cuts
            res.append(ChainMap({"dataset": newds}, spec))
        return res

    def produce_theory_database(self):
        """Produces path to the theory.db file"""
        return self.loader.theorydb_file

    def produce_combined_shift_and_theory_dataspecs(self, theoryconfig, shiftconfig):
        total_dataspecs = theoryconfig["dataspecs"] + shiftconfig["dataspecs"]
        matched_datasets = self.produce_matched_datasets_from_dataspecs(total_dataspecs)
        for ns in matched_datasets:
            ns["dataspecs"] = self.produce_dataspecs_with_matched_cuts(ns["dataspecs"])
        new_theoryconfig = []
        new_shiftconfig = []
        len_th = len(theoryconfig["dataspecs"])
        for s in matched_datasets:
            new_theoryconfig.append(ChainMap({"dataspecs": s["dataspecs"][:len_th]}, s))
            new_shiftconfig.append(ChainMap({"dataspecs": s["dataspecs"][len_th:]}, s))
        return {
            "shiftconfig": {"dataspecs": new_shiftconfig, "original": shiftconfig},
            "theoryconfig": {"dataspecs": new_theoryconfig, "original": theoryconfig},
        }

    # TODO: Worth it to do some black magic to not pass params explicitly?
    # Note that `parse_experiments` doesn't exist yet.
    def parse_reweighting_experiments(
        self, experiments, *, theoryid, use_cuts, fit=None
    ):
        """A list of experiments to be used for reweighting."""
        return self.parse_experiments(
            experiments, theoryid=theoryid, use_cuts=use_cuts, fit=fit
        )

    def parse_t0pdfset(self, name):
        """PDF set used to generate the t0 covmat."""
        return self.parse_pdf(name)

    def parse_use_t0(self, do_use_t0: bool):
        """Whether to use the t0 PDF set to generate covariance matrices."""
        return do_use_t0

    # TODO: Find a good name for this
    def produce_t0set(self, use_t0=False, t0pdfset=None):
        """Return the t0set if use_t0 is True and None otherwise. Raises an
        error if t0 is requested but no t0set is given."""
        if use_t0:
            if not t0pdfset:
                raise ConfigError("Setting use_t0 requires specifying a valid t0pdfset")
            return t0pdfset
        else:
            return None

    @element_of("posdatasets")
    def parse_posdataset(self, posset: dict, *, theoryid):
        """An observable used as positivity constrain in the fit.
        It is a mapping containing 'dataset' and 'poslambda'."""
        bad_msg = (
            "posset must be a mapping with a name ('dataset') and "
            "a float multiplier(poslambda)"
        )

        theoryno, theopath = theoryid
        try:
            name = posset["dataset"]
            poslambda = float(posset["poslambda"])
        except KeyError as e:
            raise ConfigError(bad_msg, e.args[0], posset.keys()) from e
        except ValueError as e:
            raise ConfigError(bad_msg) from e

        try:
            return self.loader.check_posset(theoryno, name, poslambda)
        except FileNotFoundError as e:
            raise ConfigError(e) from e

    def produce_posdatasets(self, positivity):
        if not isinstance(positivity, dict) or "posdatasets" not in positivity:
            raise ConfigError(
                "Failed to get 'posdatasets' from positivity. "
                "Expected that key to be present."
            )
        return positivity["posdatasets"]

    @element_of("integdatasets")
    def parse_integdataset(self, integset: dict, *, theoryid):
        """An observable corresponding to a PDF in the evolution basis,
        used as integrability constrain in the fit. It is a mapping containing 'dataset' and 'poslambda'."""
        bad_msg = (
            "integset must be a mapping with a name ('dataset') and "
            "a float multiplier(poslambda)"
        )

        theoryno, theopath = theoryid
        try:
            name = integset["dataset"]
            poslambda = float(integset["poslambda"])
        except KeyError as e:
            raise ConfigError(bad_msg, e.args[0], integset.keys()) from e
        except ValueError as e:
            raise ConfigError(bad_msg) from e
        # use the same underlying c++ code as the positivity observables
        try:
            return self.loader.check_posset(theoryno, name, poslambda)
        except FileNotFoundError as e:
            raise ConfigError(e) from e

    def produce_integdatasets(self, integrability):
        if not isinstance(integrability, dict) or "integdatasets" not in integrability:
            raise ConfigError(
                "Failed to get 'integdatasets' from integrability. "
                "Expected that key to be present."
            )
        return integrability["integdatasets"]

    def produce_reweight_all_datasets(self, experiments):
        ret = []
        for experiment in experiments:
            for dsinput, dataset in zip(experiment, experiment.datasets):
                single_exp = DataGroupSpec(
                    experiment.name, datasets=[dataset], dsinputs=[dsinput]
                )
                ret.append(
                    {"reweighting_experiments": [single_exp], "dataset_input": dsinput}
                )
        return ret

    """
    def produce_theoryid(self, theory):
        if not isinstance(theory, dict) or 'theoryid' not in theory:
            raise ConfigError("Failed to get 'theoryid' from 'theory'. "
                              "Expected that key to be present.")
        return theory['theoryid']
    """

    def produce_pdf_id(self, pdf) -> str:
        """Return a string containing the PDF's LHAPDF ID"""
        return pdf.name

    def produce_fit_id(self, fit) -> str:
        """Return a string containing the ID of the fit"""
        return fit.name

    @element_of("lumi_channels")
    def parse_lumi_channel(self, ch: str):
        if ch not in LUMI_CHANNELS:
            raise ConfigError(
                "lumi_channel not understood: %s" % ch,
                ch,
                alternatives=LUMI_CHANNELS,
                display_alternatives="all",
            )
        return ch

    def produce_all_lumi_channels(self):
        return {"lumi_channels": self.parse_lumi_channels(list(LUMI_CHANNELS))}

    @configparser.explicit_node
    def produce_nnfit_theory_covmat(
        self,
        use_thcovmat_in_sampling: bool,
        use_thcovmat_in_fitting: bool,
        thcovmat_type: str = "full",
    ):
        """
        Return the theory covariance matrix used in the fit.
        By default it is set to be the full one, the user can
        set it to be block-diagonal or diagonal, based on the
        value of ``thcovmat_type``. The possible options are:

        ``thcovmat_type = "full"`` (default):
            Include all correlations. The covarance matrix is
            computed using ``theory_covmat_custom``.

        ``thcovmat_type = "diagonal"``:
            Only diagonal entries are computes included. The
            covariance matrix is computed using
            ``theory_diagonal_covmat``.

        ``thcovmat_type = "blockdiagonal"``:
            Only correlations by process type are included.
            The covariance matrix is computed using
            ``theory_block_diag_covmat``.
        """
        valid_type = {"full", "blockdiagonal", "diagonal"}
        if thcovmat_type not in valid_type:
            raise ConfigError(
                f"Invalid thcovmat_type setting: '{valid_type}'.",
                thcovmat_type,
                valid_type,
            )

        from validphys.theorycovariance.construction import theory_covmat_custom
        from validphys.theorycovariance.construction import theory_diagonal_covmat
        from validphys.theorycovariance.construction import theory_block_diag_covmat

        if thcovmat_type == "full":
            f = theory_covmat_custom
        if thcovmat_type == "diagonal":
            f = theory_diagonal_covmat
        if thcovmat_type == "blockdiagonal":
            f = theory_block_diag_covmat

        @functools.wraps(f)
        def res(*args, **kwargs):
            return f(*args, **kwargs)

        # Set this to get the same filename regardless of the action.
        res.__name__ = "theory_covmat"
        return res

    def produce_fitthcovmat(
        self, use_thcovmat_if_present: bool = False, fit: (str, type(None)) = None
    ):
        """If a `fit` is specified and `use_thcovmat_if_present` is `True` then returns the
        corresponding covariance matrix for the given fit if it exists. If the fit doesn't have a
        theory covariance matrix then returns `False`.
        """
        if not isinstance(use_thcovmat_if_present, bool):
            raise ConfigError(
                "use_thcovmat_if_present should be a boolean, by default it is False"
            )

        if use_thcovmat_if_present and not fit:
            raise ConfigError(
                "`use_thcovmat_if_present` was true but no `fit` was specified."
            )

        if use_thcovmat_if_present and fit:
            try:
                thcovmat_present = fit.as_input()["theorycovmatconfig"][
                    "use_thcovmat_in_fitting"
                ]
            except KeyError:
                # assume covmat wasn't used and fill in key accordingly but warn user
                log.warning(
                    "use_thcovmat_if_present was true but the flag "
                    "`use_thcovmat_in_fitting` didn't exist in the runcard for "
                    f"{fit.name}. Theory covariance matrix will not be used "
                    "in any statistical estimators."
                )
                thcovmat_present = False

        if use_thcovmat_if_present and thcovmat_present:
            # Expected path of covmat hardcoded
            covmat_path = (
                fit.path
                / "tables"
                / "datacuts_theory_theorycovmatconfig_theory_covmat.csv"
            )
            if not os.path.exists(covmat_path):
                raise ConfigError(
                    "Fit appeared to use theory covmat in fit but the file was not at the "
                    f"usual location: {covmat_path}."
                )
            fit_theory_covmat = ThCovMatSpec(covmat_path)
        else:
            fit_theory_covmat = None
        return fit_theory_covmat

    def parse_speclabel(self, label: (str, type(None))):
        """A label for a dataspec. To be used in some plots"""
        return label

    @element_of("fitdeclarations")
    def parse_fitdeclaration(self, label: str):
        """Used to guess some informtion from the fit name, without having
        to download it. This is meant to be used with other providers like
        e.g.:

        {@with fits_as_from_fitdeclarations::fits_name_from_fitdeclarations@}
        {@ ...do stuff... @}
        {@endwith@}
        """
        return label

    def produce_all_commondata(self):
        """produces all commondata using the loader function """
        ds_names = self.loader.available_datasets
        ds_inputs = [self.parse_dataset_input({"dataset": ds}) for ds in ds_names]
        cd_out = [
            self.produce_commondata(dataset_input=ds_input) for ds_input in ds_inputs
        ]
        return cd_out

    def parse_groupby(self, grouping: str):
        """parses the groupby key and checks it is an allowed grouping"""
        # TODO: think if better way to do this properly
        if grouping not in ["experiment", "nnpdf31_process"]:
            raise ConfigError(
                f"Grouping not available: {grouping}, did you spell it " "correctly?"
            )
        return grouping

    def parse_norm_threshold(self, val: (numbers.Number, type(None))):
        """The threshold to use for covariance matrix normalisation, sets
        the maximum l2 norm of the inverse covariance matrix, by clipping
        smallest eigenvalues

        If norm_threshold is set to None, then no covmat regularization is
        performed

        """
        if val is not None:
            if val <= 0:
                raise ConfigError("norm_threshold must be greater than zero.")
            log.info(f"Regularizing covariance matrices with norm threshold: {val}")
        return val

    def produce_no_covmat_reg(self):
        """explicitly set norm_threshold to None so that no covariance matrix
        regularization is performed

        """
        return {"norm_threshold": None}

    @configparser.record_from_defaults
    def parse_default_filter_rules(self, spec: (str, type(None))):
        return spec

    def load_default_default_filter_rules(self, spec):
        import validphys.cuts.lockfiles

        lock_token = "_filters.lock.yaml"
        try:
            return yaml.safe_load(
                read_text(validphys.cuts.lockfiles, f"{spec}{lock_token}")
            )
        except FileNotFoundError as e:
            alternatives = [
                el.strip(lock_token)
                for el in contents(validphys.cuts.lockfiles)
                if el.endswith(lock_token)
            ]
            raise ConfigError(
                f"Default filter rules not found: {spec}",
                bad_item=spec,
                alternatives=alternatives,
                display_alternatives="all",
            )

    def parse_filter_rules(self, filter_rules: (list, type(None))):
        """A list of filter rules. See https://docs.nnpdf.science/vp/filters.html
        for details on the syntax"""
        log.warning("Overwriting filter rules")
        return filter_rules

    def parse_default_filter_rules_recorded_spec_(self, spec):
        """This function is a hacky fix for parsing the recorded spec
        of filter rules. The reason we need this function is that without
        it reportengine detects a conflict in the `dataset` key.
        """
        return spec

    def produce_rules(
        self,
        theoryid,
        use_cuts,
        defaults,
        default_filter_rules=None,
        filter_rules=None,
        default_filter_rules_recorded_spec_=None,
    ):

        """Produce filter rules based on the user defined input and defaults."""
        from validphys.filters import (
            Rule,
            RuleProcessingError,
            default_filter_rules_input,
        )

        theory_parameters = theoryid.get_description()

        if filter_rules is None:
            # Don't bother loading the rules if we are not using them.
            if use_cuts is not CutsPolicy.INTERNAL:
                return None
            if default_filter_rules_recorded_spec_ is not None:
                filter_rules = default_filter_rules_recorded_spec_[default_filter_rules]
            else:
                filter_rules = default_filter_rules_input()

        try:
            rule_list = [
                Rule(
                    initial_data=i,
                    defaults=defaults,
                    theory_parameters=theory_parameters,
                    loader=self.loader,
                )
                for i in filter_rules
            ]
        except RuleProcessingError as e:
            raise ConfigError(f"Error Processing filter rules: {e}") from e

        return rule_list

    @configparser.record_from_defaults
    def parse_default_filter_settings(self, spec: (str, type(None))):
        return spec

    def load_default_default_filter_settings(self, spec):
        import validphys.cuts.lockfiles

        lock_token = "_defaults.lock.yaml"
        try:
            return yaml.safe_load(
                read_text(validphys.cuts.lockfiles, f"{spec}{lock_token}")
            )
        except FileNotFoundError as e:
            alternatives = alternatives = [
                el.strip(lock_token)
                for el in contents(validphys.cuts.lockfiles)
                if el.endswith(lock_token)
            ]
            raise ConfigError(
                f"Default filter settings not found: {spec}",
                bad_item=spec,
                alternatives=alternatives,
                display_alternatives="all",
            )

    def parse_filter_defaults(self, filter_defaults: (dict, type(None))):
        """A mapping containing the default kinematic limits to be used when
        filtering data (when using internal cuts).
        Currently these limits are ``q2min`` and ``w2min``.
        """
        log.warning("Overwriting filter defaults")
        return filter_defaults

    def produce_defaults(
        self,
        q2min=None,
        w2min=None,
        default_filter_settings=None,
        filter_defaults={},
        default_filter_settings_recorded_spec_=None,
    ):
        """Produce default values for filters taking into account both the
        values of ``q2min`` and ` `w2min`` defined at namespace
        level and those inside a ``filter_defaults`` mapping.
        """
        from validphys.filters import default_filter_settings_input

        if (
            q2min is not None
            and "q2min" in filter_defaults
            and q2min != filter_defaults["q2min"]
        ):
            raise ConfigError("q2min defined multiple times with different values")
        if (
            w2min is not None
            and "w2min" in filter_defaults
            and w2min != filter_defaults["w2min"]
        ):
            raise ConfigError("w2min defined multiple times with different values")

        if default_filter_settings_recorded_spec_ is not None:
            filter_defaults = default_filter_settings_recorded_spec_[
                default_filter_settings
            ]
            # If we find recorded specs return immediately and don't read q2min and w2min
            # from runcard
            return filter_defaults
        elif not filter_defaults:
            filter_defaults = default_filter_settings_input()
            defaults_loaded = True
        else:
            defaults_loaded = False

        if q2min is not None and defaults_loaded:
            log.warning("Using q2min from runcard")
            filter_defaults["q2min"] = q2min

        if w2min is not None and defaults_loaded:
            log.warning("Using w2min from runcard")
            filter_defaults["w2min"] = w2min

        return filter_defaults

    def produce_data(
        self,
        data_input,
        *,
        group_name="data",
    ):
        """A set of datasets where correlated systematics are taken
        into account
        """
        datasets = []
        for dsinp in data_input:
            with self.set_context(ns=self._curr_ns.new_child({"dataset_input": dsinp})):
                datasets.append(self.parse_from_(None, "dataset", write=False)[1])

        return DataGroupSpec(name=group_name, datasets=datasets, dsinputs=data_input)

    def _parse_data_input_from_(
        self,
        parse_from_value: (str, type(None)),
        additional_context: (dict, type(None)) = None,
    ):
        """Function which parses the ``data_input`` from a namespace. Usage
        is similar to :py:meth:`self.parse_from_` except this function bridges
        the gap between the new and old way of specifying data.

        First it attempts to parse ``dataset_inputs`` from the namespace
        specified by ``parse_from_value``, for more information see
        :py:meth:`self.parse_from_`. If that fails then attempt to parse
        ``experiments``. If both should fail then raise the first exception
        encountered from the second, so that the cause can be surface in
        ``debug`` mode.

        Parameters
        ----------
        parse_from_value: str, None
            value which will be passed to :py:meth:`self.parse_from_`. If None
            then parses from the current namespace but can also be another
            input resource which can be resolved as a ``dict``.
        additional_context: dict, None
            additional context to update the namespace specified by
            ``parse_from_value``.
            In the case of this function, if ``experiments`` needs to be parsed
            then it has the additional requirements of ``theoryid`` and
            ``use_cuts`` which should either already be present in
            ``parse_from_value`` or can be passed as a ``dict`` using this
            parameter i.e ``additional_context={"theoryid": 53}``.

        """
        with self.set_context(ns=self._curr_ns.new_child(additional_context)):
                # new fits have dataset_inputs, old fits have experiments
            data_key = "dataset_inputs"
            try:
                _, data_val = self.parse_from_(parse_from_value, data_key, write=False)
            except ConfigError as e:
                data_key = "experiments"
                log.warning(
                    "`experiments` has been deprecated, specify data using `dataset_inputs`. "
                    "Any grouping defined by `experiments` is being ignored."
                )
                # We need to make theoryid available if using experiments
                try:
                    _, experiments = self.parse_from_(parse_from_value, data_key, write=False)
                    data_val = NSList([
                        dsinput for experiment in experiments for dsinput in experiment.dsinputs
                    ], nskey='dataset_input')
                except ConfigError as inner_error:
                    log.error(inner_error)
                    raise e from inner_error
        return data_val

    def produce_data_input(self):
        """Produce the ``data_input`` which is a flat list of ``dataset_input`` s.
        This production rule handles the backwards compatibility with old datasets
        which specify ``experiments`` in the runcard.

        """
        # parse from current namespace with no additional context.
        return self._parse_data_input_from_(None)

    def parse_metadata_group(self, group: str):
        """User specified key to group data by. The key must exist in the
        PLOTTING file for example `experiment`
        """
        return group

    @record_from_defaults
    def parse_data_grouping(self, key):
        """a key which indicates which default grouping to use. Mainly for
        internal use. It allows the default grouping of experiment to be applied
        to runcards which don't specify `metadata_group` without there being
        a namespace conflict in the lockfile

        """
        return key

    def load_default_data_grouping(self, spec):
        """Load the default grouping of data"""
        # slightly superfluous, only one default at present but perhaps
        # somebody will want to add to this at some point e.g for th. uncertainties
        allowed = {
            "standard_report": "experiment",
        }
        return allowed[spec]

    def produce_processed_data_grouping(
        self, data_grouping=None, data_grouping_recorded_spec_=None
    ):
        """Process the data_grouping key from the runcard, or lockfile. If
        `data_grouping_recorded_spec_` is present then its value is taken, and
        the runcard is assumed to be a lockfile.

        If data_grouping is None, then fall back to old behaviour of grouping
        by experiment.

        Else, the user can specfiy their own grouping, for example metadata_process.
        """
        if data_grouping is None:
            # fallback to old default behaviour, but still record to lockfile
            data_grouping = self.parse_data_grouping("standard_report")
        if data_grouping_recorded_spec_ is not None:
            return data_grouping_recorded_spec_[data_grouping]
        return self.load_default_data_grouping(data_grouping)

    def produce_processed_metadata_group(
        self, processed_data_grouping, metadata_group=None
    ):
        """Expose the final data grouping result. Either metadata_group is
        specified by user, in which case uses `processed_data_grouping` which
        is experiment by default.
        """
        if metadata_group is None:
            return processed_data_grouping
        return metadata_group


    def produce_group_dataset_inputs_by_metadata(
        self, data_input, processed_metadata_group,
    ):
        """Take the data and the processed_metadata_group key and attempt
        to group the data, returns a list where each element specifies the data_input
        for a single group and the group_name
        """
        res = defaultdict(list)
        for dsinput in data_input:
            cd = self.produce_commondata(dataset_input=dsinput)
            try:
                res[getattr(get_info(cd), processed_metadata_group)].append(dsinput)
            except AttributeError as e:
                raise ConfigError(
                    f"Unable to find key: {processed_metadata_group} in {cd.name} "
                    "PLOTTING file.",
                    bad_item=processed_metadata_group,
                    alternatives=get_info(cd).__dict__,
                ) from e
        return [
            {"data_input": NSList(group, nskey="dataset_input"), "group_name": name}
            for name, group in res.items()
        ]


    def produce_group_dataset_inputs_by_experiment(self, data_input):
        return self.produce_group_dataset_inputs_by_metadata(data_input, "experiment")

    def produce_scale_variation_theories(self, theoryid, point_prescription):
        """Produces a list of theoryids given a theoryid at central scales and a point
           prescription. The options for the latter are '3 point', '5 point', '5bar point', '7 point'
           and '9 point'. Note that these are defined in arXiv:1906.10698. This hard codes the
           theories needed for each prescription to avoid user error."""
        pp = point_prescription
        th = theoryid.id

        lsv = yaml.safe_load(
            read_text(validphys.scalevariations, "scalevariationtheoryids.yaml")
        )

        scalevarsfor_list = lsv["scale_variations_for"]
        # Allowed central theoryids
        cent_thids = [
            str(scalevarsfor_dict["theoryid"])
            for scalevarsfor_dict in scalevarsfor_list
        ]

        if th not in cent_thids:
            valid_thids = ", ".join(cent_thids)
            raise ConfigError(
                "Scale variations are not currently defined for this central theoryid. It is "
                + f"currently only possible to use one of the following as the central theory: {valid_thids}. "
                + "Please use one of these instead if you wish to include theory uncertainties here."
            )

        # Find scales that correspond to this point prescription
        pp_scales_dict = yaml.safe_load(
            read_text(validphys.scalevariations, "pointprescriptions.yaml")
        )

        try:
            scales = pp_scales_dict[pp]
        except KeyError:
            valid_pps = ", ".join(pp_scales_dict.keys())
            raise ConfigError(
                "Scale variations are not currently defined for this point prescription. This "
                + "configuration only works when 'point_prescription' is equal to one of the "
                + f"following: {valid_pps}. Please use one of these instead if you wish to "
                + "include theory uncertainties here."
            )

        # Get dictionary containing theoryid and variations for central theory from runcard
        for scalevarsfor_dict in scalevarsfor_list:
            if scalevarsfor_dict["theoryid"] == int(th):
                theoryid_variations = scalevarsfor_dict

        # Find theoryids for given point prescription for given central theoryid
        try:
            thids = [theoryid_variations["variations"][scale] for scale in scales]
        except KeyError:
            available_scales = list(theoryid_variations["variations"])
            missing_scales = []
            for scale in scales:
                if scale not in available_scales:
                    missing_scales.append(scale)
            missing_scales_string = ", ".join(missing_scales)
            raise ConfigError(
                "For this central theoryid, the requested point prescription is not currently "
                + "available. To use this point prescription for this central theoryid, theoryids "
                + "that correspond to the following scale choices must be created and added to "
                + "validphys2/src/validphys/scalevariations/scalevariationtheoryids.yaml: "
                + f"(k_F, k_R) = {missing_scales_string}."
            )

        # Check each theory is loaded
        theoryids = [self.loader.check_theoryID(thid) for thid in thids]
        # NSList needs to be used for theoryids to be recognised as a namespace
        return {"theoryids": NSList(theoryids, nskey="theoryid")}


    @configparser.explicit_node
    def produce_filter_data(self, fakedata: bool = False, theorycovmatconfig=None):
        """Set the action used to filter the data to filter either real or
        closure data. If the closure data filter is being used and if the
        theory covariance matrix is not being closure tested then filter
        data by experiment for efficiency"""
        import validphys.filters
        if not fakedata:
            return validphys.filters.filter_real_data
        else:
            if (
                theorycovmatconfig is not None and
                theorycovmatconfig.get("use_thcovmat_in_sampling")
            ):
                # NOTE: By the time we run theory covmat closure tests,
                # hopefully the generation of pseudodata will be done in python.
                raise ConfigError(
                    "Generating closure test data which samples from the theory "
                    "covariance matrix has not been implemented yet."
                )
            return validphys.filters.filter_closure_data_by_experiment

    @configparser.explicit_node
    def produce_total_chi2_data(self, fitthcovmat):
        """If there is no theory covmat for the fit, then calculate the
        total chi2 by summing the chi2 from each experiment.
        """
        import validphys.results

        if fitthcovmat is None:
            return validphys.results.total_chi2_data_from_experiments
        return validphys.results.dataset_inputs_abs_chi2_data

    @configparser.explicit_node
    def produce_total_phi_data(self, fitthcovmat):
        """If there is no theory covmat for the fit, then calculate the total
        phi using contributions from each experiment.
        """
        import validphys.results

        if fitthcovmat is None:
            return validphys.results.total_phi_data_from_experiments
        return validphys.results.dataset_inputs_phi_data



class Config(report.Config, CoreConfig, ParamfitsConfig):
    """The effective configuration parser class."""

    pass
