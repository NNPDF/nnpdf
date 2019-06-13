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

from collections import ChainMap
from collections.abc import Mapping, Sequence

from reportengine import configparser
from reportengine.environment import Environment, EnvironmentError_
from reportengine.configparser import ConfigError, element_of, _parse_func
from reportengine.helputils import get_parser_type
from reportengine import report

from validphys.core import (ExperimentSpec, DataSetInput, ExperimentInput,
                            CutsPolicy, MatchedCuts, ThCovMatSpec)
from validphys.loader import (Loader, LoaderError ,LoadFailedError, DataNotFoundError,
                              PDFNotFound, FallbackLoader)
from validphys.gridvalues import LUMI_CHANNELS

from validphys.paramfits.config import ParamfitsConfig

from validphys.theorycovariance.theorycovarianceutils import process_lookup
from validphys.plotoptions import get_info

log = logging.getLogger(__name__)

class Environment(Environment):
    """Container for information to be filled at run time"""

    def __init__(self,*, datapath=None, resultspath=None, this_folder=None, net=True,
                 upload=False, **kwargs):
        if this_folder:
            self.this_folder = pathlib.Path(this_folder)

        if net:
            loader_class = FallbackLoader
        else:
            loader_class = Loader

        try:
            self.loader = loader_class(datapath, resultspath)
        except LoaderError as e:
            log.error("Failed to find the paths. These are configured "
                      "in the nnprofile settings or with the "
                      "--datapath and --resultpath options.")
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
        keydiff =  item.keys() - {'id', 'label'}

        if  keydiff or not 'id' in item:
            unrecognized = f' Unrecognized: {keydiff}' if keydiff else ''
            raise ConfigError(f"'{item}' must be a single id, or a mapping "
                              f"with keys 'id', 'label.{unrecognized}'")
        id = item['id']
        val = f(self, id, **kwargs)
        if 'label' in item:
            val.label = str(item['label'])
        return val

    currsig = inspect.signature(parse_func)
    origsig = inspect.signature(f)
    parse_func = functools.wraps(f)(parse_func)

    params = [*list(currsig.parameters.values())[:2],
              *list(origsig.parameters.values())[2:]]

    parse_func.__signature__ = inspect.Signature(
                 parameters=params)

    labeldoc =  ((" Either just an id %s, or a mapping "
                          "with 'id' and 'label'.") %
                          (get_parser_type(f),))
    if parse_func.__doc__ is None:
        parse_func.__doc__ = labeldoc
    else:
        parse_func.__doc__ += labeldoc

    return parse_func

class CoreConfig(configparser.Config):


    @property
    def loader(self):
        return self.environment.loader

    @element_of('pdfs')
    @_id_with_label
    def parse_pdf(self, name:str):
        """A PDF set installed in LHAPDF."""
        try:
            pdf = self.loader.check_pdf(name)
        except PDFNotFound as e:
            raise ConfigError("Bad PDF: {} not installed".format(name), name,
                          self.loader.available_pdfs) from e
        except LoaderError as e:
            raise ConfigError(e) from e

        #Check that we know how to compute errors
        try:
            pdf.nnpdf_error
        except NotImplementedError as e:
            raise ConfigError(str(e))
        return pdf


    @element_of('theoryids')
    @_id_with_label
    def parse_theoryid(self, theoryID: (str, int)):
        """A number corresponding to the database theory ID where the
        corresponding theory folder is installed in te data directory."""
        try:
            return self.loader.check_theoryID(theoryID)
        except LoaderError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    def parse_use_cuts(self, use_cuts: (bool, str)):
        """Whether to filter the points based on the cuts applied in the fit,
        or the whole data in the dataset. The possible options are:

        - internal: Calculate the cuts based on the existing rules. This is
          the default.

        - fromfit: Read the cuts stored in the fit.

        - nocuts: Use the whole dataset.
        """
        #The lower is an aesthetic preference...
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
                f"equivalent to '{res}'.")
        elif isinstance(use_cuts, str) and use_cuts in valid_cuts:
            res = CutsPolicy(use_cuts)
        else:
            raise ConfigError(f"Invalid use_cuts setting: '{use_cuts}'.",
                              use_cuts, valid_cuts)

        return res


    #TODO: load fit config from here
    @element_of('fits')
    @_id_with_label
    def parse_fit(self, fit:str):
        """A fit in the results folder, containing at least a valid filter result."""
        try:
            return self.loader.check_fit(fit)
        except LoadFailedError as e:
            raise ConfigError(str(e), fit ,self.loader.available_fits)

    def produce_fitcontext(self, fit):
        """Set PDF, theory ID and experiments from the fit config"""

        _, pdf         = self.parse_from_('fit', 'pdf', write=False)
        _, theory      = self.parse_from_('fit', 'theory', write=False)

        #TODO: parse we need multilevel from to do theoruid nicely
        thid = theory['theoryid']

        #We need to make theoryid available to parse the experiments
        with self.set_context(ns=self._curr_ns.new_child({'theoryid':thid})):
            _, experiments = self.parse_from_('fit', 'experiments', write=False)



        return {'pdf': pdf, 'theoryid':thid, 'experiments': experiments}

    def produce_fitinputcontext(self, fit):
        """Like ``fitcontext`` but without setting the PDF"""


        _, theory      = self.parse_from_('fit', 'theory', write=False)
        thid = theory['theoryid']

        #We need to make theoryid available to parse the experiments
        with self.set_context(ns=self._curr_ns.new_child({'theoryid':thid})):
            _, experiments = self.parse_from_('fit', 'experiments', write=False)

        return {'theoryid':thid, 'experiments': experiments}

    def produce_fitpdf(self, fit):
        """Like ``fitcontext`` only setting the PDF"""
        with self.set_context(ns=self._curr_ns.new_child({'fit':fit})):
            _, pdf = self.parse_from_('fit', 'pdf', write=False)
        return {'pdf': pdf}

    def produce_fitunderlyinglaw(self, fit):
        """Reads closuretest: fakepdf from fit config file and passes as
        pdf
        """
        with self.set_context(ns=self._curr_ns.new_child({'fit':fit})):
            _, datacuts = self.parse_from_('fit', 'closuretest', write=False)
        underlyinglaw = self.parse_pdf(datacuts['fakepdf'])
        return {'pdf': underlyinglaw}

    def produce_fitpdfandbasis(self, fit):
        """ Set the PDF and basis from the fit config. """
        with self.set_context(ns=self._curr_ns.new_child({'fit':fit})):
            _, pdf = self.parse_from_('fit', 'pdf', write=False)
            _, fitting = self.parse_from_('fit', 'fitting', write=False)
        basis = fitting['fitbasis']
        return {'pdf': pdf, 'basis':basis}

    @element_of('dataset_inputs')
    def parse_dataset_input(self, dataset:Mapping):
        """The mapping that corresponds to the dataset specifications in the
        fit files"""
        known_keys = {'dataset', 'sys', 'cfac', 'frac', 'weight'}
        try:
            name = dataset['dataset']
            if not isinstance(name, str):
                raise ConfigError(f"'dataset' must be a string, not {type(name)}")
        except KeyError:
            raise ConfigError("'dataset' must be a mapping with "
                              "'dataset' and 'sysnum'")


        sysnum = dataset.get('sys')
        cfac = dataset.get('cfac', tuple())
        frac = dataset.get('frac', 1)
        if  not isinstance(frac, numbers.Real):
            raise ConfigError(f"'frac' must be a number, not '{frac}'")
        if frac < 0 or frac > 1:
            raise ConfigError(f"'frac' must be between 0 and 1 not '{frac}'")
        weight = dataset.get('weight', 1)
        if  not isinstance(weight, numbers.Real):
            raise ConfigError(f"'weight' must be a number, not '{weight}'")
        if weight < 0:
            raise ConfigError(f"'weight' must be greater than zero not '{weight}'")
        kdiff = dataset.keys() - known_keys
        for k in kdiff:
            #Abuse ConfigError to get the suggestions.
            log.warninig(ConfigError(f"Key '{k}' in dataset_input not known.", k, known_keys))
        return DataSetInput(name=name, sys=sysnum, cfac=cfac, frac=frac,
                weight=weight)

    def parse_use_fitcommondata(self, do_use: bool):
        """Use the commondata files in the fit instead of those in the data
        directory."""
        return do_use

    def produce_commondata(self,
                           *,
                           dataset_input,
                           use_fitcommondata=False,
                           fit=None):
        """Produce a CommondataSpec from a dataset input"""

        name = dataset_input.name
        sysnum = dataset_input.sys
        try:
            return self.loader.check_commondata(
                setname=name,
                sysnum=sysnum,
                use_fitcommondata=use_fitcommondata,
                fit=fit)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name,
                              self.loader.available_datasets) from e
        except LoadFailedError as e:
            raise ConfigError(e) from e

    def produce_cuts(self,
                     *,
                     commondata,
                     use_cuts,
                     fit=None,
                     q2min: (numbers.Number, type(None)) = None,
                     w2min: (numbers.Number, type(None)) = None,
                     theoryid=None):
        """Obtain cuts for a given dataset input, based on the
        appropriate policy."""
        #TODO: Put this bit of logic into loader.check_cuts
        if use_cuts is CutsPolicy.NOCUTS:
            return None
        elif use_cuts is CutsPolicy.FROMFIT:
            if not fit:
                raise ConfigError(
                    "Setting 'use_cuts' to 'fromfit' requires "
                    "specifying a fit on which filter "
                    "has been executed, e.g.\nfit : NNPDF30_nlo_as_0118")
            name = commondata.name
            try:
                return self.loader.check_fit_cuts(name, fit)
            except LoadFailedError as e:
                raise ConfigError(e) from e
        elif use_cuts is CutsPolicy.INTERNAL:
            if not isinstance(q2min, numbers.Number):
                raise ConfigError("q2min must be specified for internal cuts")
            if not isinstance(w2min, numbers.Number):
                raise ConfigError("w2min must be specified for internal cuts")
            if not theoryid:
                raise ConfigError(
                    "thoeryid must be specified for internal cuts")
            return self.loader.check_internal_cuts(commondata, theoryid, q2min, w2min)
        elif use_cuts is CutsPolicy.FROM_CUT_INTERSECTION_NAMESPACE:
            cut_list = []
            _, nss = self.parse_from_(None, 'cuts_intersection_spec', write=False)
            self._check_dataspecs_type(nss)
            if not nss:
                raise ConfigError("'cuts_intersection_spec' must contain at least one namespace.")

            ns_cut_inputs = {'commondata': commondata,
                             'use_cuts': CutsPolicy.INTERNAL}
            for ns in nss:
                with self.set_context(ns=self._curr_ns.new_child({**ns, **ns_cut_inputs})):
                    _, nscuts = self.parse_from_(None, 'cuts', write=False)
                    cut_list.append(nscuts)
            ndata = commondata.ndata
            return MatchedCuts(cut_list, ndata=ndata)

        raise TypeError("Wrong use_cuts")


    def produce_dataset(self,
                        *,
                        dataset_input,
                        theoryid,
                        cuts,
                        use_fitcommondata=False,
                        fit=None,
                        check_plotting: bool = False):
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
                name=name,
                sysnum=sysnum,
                theoryid=theoryid,
                cfac=cfac,
                cuts=cuts,
                frac=frac,
                use_fitcommondata=use_fitcommondata,
                fit=fit,
                weight=weight)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets)

        except LoadFailedError as e:
            raise ConfigError(e)

        if check_plotting:
            from validphys.plotoptions import get_info
            #normalize=True should check for more stuff
            get_info(ds, normalize=True)
            if not ds.commondata.plotfiles:
                log.warning(f"Plotting files not found for: {ds}")
        return ds

    @configparser.element_of('experiments')
    def parse_experiment(
            self,
            experiment: dict,
            *,
            theoryid,
            use_cuts,
            fit=None,
            check_plotting: bool = False,
            use_fitcommondata=False,
            q2min: (numbers.Number, type(None)) = None,
            w2min: (numbers.Number, type(None)) = None,
    ):
        """A set of datasets where correlated systematics are taken
           into account. It is a mapping where the keys are the experiment
           name 'experiment' and a list of datasets."""
        try:
            name, datasets = experiment['experiment'], experiment['datasets']
        except KeyError as e:
            raise ConfigError(
                "'experiment' must be a mapping with "
                "'experiment' and 'datasets', but %s is missing" % e)

        dsinputs = [self.parse_dataset_input(ds) for ds in datasets]
        cds = [self.produce_commondata(
                dataset_input=dsinp,
                use_fitcommondata=use_fitcommondata,
                fit=fit) for dsinp in dsinputs]
        cutinps = [
            self.produce_cuts(
                commondata=cd,
                use_cuts=use_cuts,
                q2min=q2min,
                w2min=w2min,
                fit=fit,
                theoryid=theoryid,
                ) for cd in cds
        ]

        #autogenerated func, from elemet_of
        datasets = [
            self.produce_dataset(
                dataset_input=dsinp,
                theoryid=theoryid,
                cuts=cuts,
                fit=fit,
                check_plotting=check_plotting,
                use_fitcommondata=use_fitcommondata)
            for (dsinp, cuts) in zip(dsinputs, cutinps)
        ]

        return ExperimentSpec(name=name, datasets=datasets, dsinputs=dsinputs)


    @configparser.element_of('experiment_inputs')
    def parse_experiment_input(self, ei:dict):
        """The mapping that corresponds to the experiment specification in the
        fit config files. Currently, this needs to be combined with
        ``experiment_from_input`` to yield an useful result."""
        try:
            name = ei['experiment']
        except KeyError as e:
            raise ConfigError(f"experiment_input must have an 'experiment' key") from e

        try:
            datasets = ei['datasets']
        except KeyError as e:
            raise ConfigError(f"experiment_input must have an 'datasets' key") from e

        return ExperimentInput(name=name, datasets=datasets)

    #TODO: Do away with the mapping and make the conversion implicitly
    def produce_experiment_from_input(self, experiment_input,theoryid, use_cuts,
            fit=None):
        """Return a mapping containing a single experiment from an experiment
        input. NOTE: This might be deprecated in the future."""
        return {'experiment':self.parse_experiment(experiment_input.as_dict(),
                theoryid=theoryid, use_cuts=use_cuts, fit=fit)}

    #TODO: Do this better and elsewhere
    @staticmethod
    def _check_dataspecs_type(dataspecs):
        if not isinstance(dataspecs, Sequence):
            raise ConfigError(
                "dataspecs should be a sequence of mappings, not "
                f"{type(dataspecs).__name__}")

        for spec in dataspecs:
            if not isinstance(spec, Mapping):
                raise ConfigError(
                    "dataspecs should be a sequence of mappings, "
                    f" but {spec} is {type(spec).__name__}")

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
                _, experiments = self.parse_from_(
                    None, 'experiments', write=False)
                names = {(e.name, ds.name, process_lookup(ds.name)): (ds, dsin)
                         for e in experiments
                         for ds, dsin in zip(e.datasets, e)}
                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))

        res = []
        for k in used_set:
            inres = {'experiment_name': k[0], 'dataset_name': k[1], 'process': k[2]}
            #TODO: Should this have the same name?
            inner_spec_list = inres['dataspecs'] = []
            for ispec, spec in enumerate(dataspecs):
                #Passing spec by referene
                d = ChainMap({
                    'dataset': all_names[ispec][k][0],
                    'dataset_input': all_names[ispec][k][1],
                }, spec)
                inner_spec_list.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x['process'], x['experiment_name'], x['dataset_name']))
        return res

    def produce_matched_positivity_from_dataspecs(self, dataspecs):
        """Like produce_matched_datasets_from_dataspecs but for positivity datasets."""
        self._check_dataspecs_type(dataspecs)
        all_names = []
        for spec in dataspecs:
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, pos = self.parse_from_(None, 'posdatasets', write=False)
                names = {(p.name):(p) for p in pos}
                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))

        res = []
        for k in used_set:
            inres = {'posdataset_name': k}
            #TODO: Should this have the same name?
            l = inres['dataspecs'] = []
            for ispec, spec in enumerate(dataspecs):
                #Passing spec by referene
                d = ChainMap({
                    'posdataset':       all_names[ispec][k],

                    },
                    spec)
                l.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x['posdataset_name']))
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
        #Can now assume we have at least one element
        cutlist = []
        dslist = []
        names = set()
        for spec in dataspecs:
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, ds = self.parse_from_(None, 'dataset', write=False)
            dslist.append(ds)
            cutlist.append(ds.cuts)
            names.add(ds.name)

        lnames = len(names)
        if lnames != 1:
            raise ConfigError("Each dataspec must have a dataset with the same"
                f"name, but got {lnames} different ones: {names}")

        ndata = ds.commondata.ndata
        matched_cuts = MatchedCuts(cutlist, ndata=ndata)
        res = []
        for spec, ds in zip(dataspecs, dslist):
            newds = copy.copy(ds)
            newds.cuts = matched_cuts
            res.append(ChainMap({'dataset': newds}, spec))
        return res

    def produce_combined_shift_and_theory_dataspecs(self, theoryconfig, shiftconfig):
        total_dataspecs = theoryconfig["dataspecs"] + shiftconfig["dataspecs"]
        matched_datasets = self.produce_matched_datasets_from_dataspecs(total_dataspecs)
        for ns in matched_datasets:
            ns["dataspecs"] = self.produce_dataspecs_with_matched_cuts(ns["dataspecs"])
        new_theoryconfig = []
        new_shiftconfig = []
        len_th = len(theoryconfig['dataspecs'])
        for s in matched_datasets:
            new_theoryconfig.append(ChainMap({"dataspecs": s['dataspecs'][:len_th]}, s))
            new_shiftconfig.append(ChainMap({"dataspecs": s['dataspecs'][len_th:]}, s))
        return {
            "shiftconfig": {"dataspecs": new_shiftconfig, "original": shiftconfig},
            "theoryconfig": {"dataspecs": new_theoryconfig, "original": theoryconfig}
        }


    #TODO: Worth it to do some black magic to not pass params explicitly?
    #Note that `parse_experiments` doesn't exist yet.
    def parse_reweighting_experiments(self, experiments, *, theoryid,
                                      use_cuts, fit=None):
        """A list of experiments to be used for reweighting."""
        return self.parse_experiments(experiments,
                                     theoryid=theoryid,
                                     use_cuts=use_cuts, fit=fit)
    def parse_t0pdfset(self, name):
        """PDF set used to generate the t0 covmat."""
        return self.parse_pdf(name)

    def parse_use_t0(self, do_use_t0:bool):
        """Whether to use the t0 PDF set to generate covariance matrices."""
        return do_use_t0

    #TODO: Find a good name for this
    def produce_t0set(self, use_t0=False, t0pdfset=None):
        """Return the t0set if use_t0 is True and None otherwise. Raises an
        error if t0 is requested but no t0set is given."""
        if use_t0:
            if not t0pdfset:
               raise ConfigError("Setting use_t0 requires specifying a valid t0pdfset")
            return t0pdfset
        else:
            return None

    @element_of('posdatasets')
    def parse_posdataset(self, posset:dict, * ,theoryid):
        """An observable used as positivity constrain in the fit.
        It is a mapping containing 'dataset' and 'poslambda'."""
        bad_msg = ("posset must be a mapping with a name ('dataset') and "
                   "a float multiplier(poslambda)")

        theoryno, theopath = theoryid
        try:
            name = posset['dataset']
            poslambda = float(posset['poslambda'])
        except KeyError as e:
            raise ConfigError(bad_msg, e.args[0], posset.keys()) from e
        except ValueError as e:
            raise ConfigError(bad_msg) from e

        try:
            return self.loader.check_posset(theoryno, name, poslambda)
        except FileNotFoundError as e:
            raise ConfigError(e) from e

    def produce_posdatasets(self, positivity):
        if not isinstance(positivity, dict) or 'posdatasets' not in positivity:
            raise ConfigError("Failed to get 'posdatasets' from positivity. "
                              "Expected that key to be present.")
        return positivity['posdatasets']

    def produce_reweight_all_datasets(self, experiments):
        ret = []
        for experiment in experiments:
            for dsinput, dataset in zip(experiment, experiment.datasets):
                single_exp = ExperimentSpec(experiment.name,
                                            datasets=[dataset],
                                            dsinputs=[dsinput])
                ret.append({'reweighting_experiments': [single_exp],
                            'dataset_input':dsinput})
        return ret

    """
    def produce_theoryid(self, theory):
        if not isinstance(theory, dict) or 'theoryid' not in theory:
            raise ConfigError("Failed to get 'theoryid' from 'theory'. "
                              "Expected that key to be present.")
        return theory['theoryid']
    """

    def produce_pdf_id(self, pdf)-> str:
        """Return a string containing the PDF's LHAPDF ID"""
        return pdf.name


    @element_of('lumi_channels')
    def parse_lumi_channel(self, ch:str):
        if ch not in LUMI_CHANNELS:
            raise ConfigError('lumi_channel not understood: %s' % ch,
                              ch, alternatives=LUMI_CHANNELS,
                              display_alternatives='all')
        return ch

    def produce_all_lumi_channels(self):
        return {'lumi_channels': self.parse_lumi_channels(list(LUMI_CHANNELS))}

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
            self, use_thcovmat_if_present: bool = False, fit: (str, type(None)) = None):
        """If a `fit` is specified and `use_thcovmat_if_present` is `True` then returns the
        corresponding covariance matrix for the given fit if it exists. If the fit doesn't have a
        theory covariance matrix then returns `False`.
        """
        if not isinstance(use_thcovmat_if_present, bool):
            raise ConfigError("use_thcovmat_if_present should be a boolean, by default it is False")

        if use_thcovmat_if_present and not fit:
            raise ConfigError("`use_thcovmat_if_present` was true but no `fit` was specified.")

        if use_thcovmat_if_present and fit:
            try:
                thcovmat_present = fit.as_input()[
                    'theorycovmatconfig']['use_thcovmat_in_fitting']
            except KeyError:
                #assume covmat wasn't used and fill in key accordingly but warn user
                log.warning("use_thcovmat_if_present was true but the flag "
                            "`use_thcovmat_in_fitting` didn't exist in the runcard for "
                            f"{fit.name}. Theory covariance matrix will not be used "
                            "in any statistical estimators.")
                thcovmat_present = False


        if use_thcovmat_if_present and thcovmat_present:
            # Expected path of covmat hardcoded
            covmat_path = (
                fit.path/'tables'/'datacuts_theory_theorycovmatconfig_theory_covmat.csv')
            if not os.path.exists(covmat_path):
                raise ConfigError(
                    "Fit appeared to use theory covmat in fit but the file was not at the "
                    f"usual location: {covmat_path}.")
            fit_theory_covmat = ThCovMatSpec(covmat_path)
        else:
            fit_theory_covmat = None
        return fit_theory_covmat

    def parse_speclabel(self, label:(str, type(None))):
        """A label for a dataspec. To be used in some plots"""
        return label

    @element_of('fitdeclarations')
    def parse_fitdeclaration(self, label:str):
        """Used to guess some informtion from the fit name, without having
        to download it. This is meant to be used with other providers like
        e.g.:

        {@with fits_as_from_fitdeclarations::fits_name_from_fitdeclarations@}
        {@ ...do stuff... @}
        {@endwith@}
        """
        return label

    def produce_fit_data_groupby_experiment(self, fit):
        """Used to produce data from the fit grouped into experiments,
        where each experiment is a group of datasets according to the experiment
        key in the plotting info file.
        """
        #TODO: consider this an implimentation detail
        from reportengine.namespaces import NSList

        with self.set_context(ns=self._curr_ns.new_child({'fit':fit})):
            _, experiments = self.parse_from_('fit', 'experiments', write=False)

        flat = (ds for exp in experiments for ds in exp.datasets)
        metaexps = [get_info(ds).experiment for ds in flat]
        res = {}
        for exp in experiments:
            for dsinput, ds in zip(exp.dsinputs, exp.datasets):
                metaexp = get_info(ds).experiment
                if metaexp in res:
                    res[metaexp].append(ds)
                else:
                    res[metaexp] = [ds]
        exps = []
        for exp in res:
            exps.append(ExperimentSpec(exp, res[exp]))

        experiments = NSList(exps, nskey='experiment')
        return {'experiments': experiments}

    def produce_fit_context_groupby_experiment(self, fit):
        """produces experiments similarly to `fit_data_groupby_experiment`
        but also sets fitcontext (pdf and theoryid)
        """
        _, pdf         = self.parse_from_('fit', 'pdf', write=False)
        _, theory      = self.parse_from_('fit', 'theory', write=False)
        thid = theory['theoryid']
        with self.set_context(ns=self._curr_ns.new_child({'theoryid':thid})):
            experiments = self.produce_fit_data_groupby_experiment(
                fit)['experiments']
        return {'pdf': pdf, 'theoryid':thid, 'experiments': experiments}



class Config(report.Config, CoreConfig, ParamfitsConfig):
    """The effective configuration parser class."""
    pass

