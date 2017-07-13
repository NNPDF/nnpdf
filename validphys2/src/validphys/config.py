# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:43:10 2016

@author: Zahari Kassabov
"""
import logging
import pathlib
import functools
import inspect
import re
from collections import Sequence, Mapping, ChainMap

from  reportengine import configparser
from reportengine.environment import Environment, EnvironmentError_
from reportengine.configparser import ConfigError, element_of, _parse_func
from reportengine.helputils import get_parser_type
from reportengine import report

from validphys.core import ExperimentSpec, DataSetInput, ExperimentInput
from validphys.loader import (Loader, LoaderError ,LoadFailedError, DataNotFoundError,
                              PDFNotFound, FallbackLoader)
from validphys import tableloader
from validphys.gridvalues import LUMI_CHANNELS

log = logging.getLogger(__name__)



class Environment(Environment):
    """Container for information to be filled at run time"""

    def __init__(self,*, datapath=None, resultspath=None, this_folder, net=True,
                 upload=False, **kwargs):

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

class Config(report.Config):


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
        except FileNotFoundError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    def parse_use_cuts(self, use_cuts:bool, *, fit=None):
        """Whether to use the filtered points in the fit, or the whole
        data in the dataset."""
        if use_cuts and not fit:
            raise ConfigError("Setting 'use_cuts' true requires "
            "specifying a fit on which filter "
            "has been executed, e.g.\nfit : NNPDF30_nlo_as_0118")
        return use_cuts

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




    @element_of('dataset_inputs')
    def parse_dataset_input(self, dataset):
        """The mapping that corresponds to the dataset specifications in the
        fit files"""
        try:
            name = dataset['dataset']
            if not isinstance(name, str):
                raise Config("'name' must be a string, not %s" % type(name))
        except KeyError:
            raise ConfigError("'dataset' must be a mapping with "
                              "'dataset' and 'sysnum'")


        sysnum = dataset.get('sys')
        cfac = dataset.get('cfac', tuple())
        return DataSetInput(name=name, sys=sysnum, cfac=cfac)

    def produce_commondata(self, *, dataset_input):
        """Produce a CommondataSpec from a dataset input"""

        name = dataset_input.name
        sysnum = dataset_input.sys
        try:
            return self.loader.check_commondata(setname=name, sysnum=sysnum)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets) from e
        except LoadFailedError as e:
            raise ConfigError(e) from e

    #TODO: Decide how should this interact with use_cuts
    def produce_cuts(self, *, dataset_input, use_cuts, fit=None):
        """Obtain cuts from a fit, for a given dataset input"""
        if not use_cuts:
            return None
        name = dataset_input.name
        try:
            return self.loader.check_cuts(name, fit)
        except LoadFailedError as e:
            raise ConfigError(e) from e


    def produce_dataset(self, *, dataset_input ,theoryid, use_cuts, fit=None,
                      check_plotting:bool=False):
        """Dataset specification from the theory and CommonData.
           Use the cuts from the fit, if provided. If check_plotting is set to
           True, attempt to lod and check the PLOTTING files
           (note this may cause a noticeable slowdown in general)."""
        name = dataset_input.name
        sysnum = dataset_input.sys
        cfac = dataset_input.cfac

        try:
            ds =  self.loader.check_dataset(name=name, sysnum=sysnum,
                                             theoryid=theoryid, cfac=cfac,
                                             use_cuts=use_cuts, fit=fit)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets)

        except LoadFailedError as e:
            raise ConfigError(e)

        if check_plotting:
            from validphys.plotoptions import get_info
            #normalize=True should check for more stuff
            get_info(ds, normalize=True)
            if not ds.commondata.plotfiles:
                log.warn("Plotting files not found for: %s" % (ds,))
        return ds


    @configparser.element_of('experiments')
    def parse_experiment(self, experiment:dict, *, theoryid, use_cuts,
                         fit=None, check_plotting:bool=False):
        """A set of datasets where correlated systematics are taken
           into account. It is a mapping where the keys are the experiment
           name 'experiment' and a list of datasets."""
        try:
            name, datasets = experiment['experiment'], experiment['datasets']
        except KeyError as e:
            raise ConfigError("'experiment' must be a mapping with "
                              "'experiment' and 'datasets', but %s is missing" % e)

        dsinputs = [self.parse_dataset_input(ds) for ds in datasets]
        #autogenerated func, from elemet_of
        datasets = [self.produce_dataset(dataset_input=dsinp, theoryid=theoryid,
                                       use_cuts=use_cuts, fit=fit,
                                       check_plotting=check_plotting)
                                       for dsinp in dsinputs]

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


    #TODO: Get rid of this
    def produce_fits_pdf_config(self, fits):
        """DO NOT USE. For internal use only,"""
        return [self.produce_fitpdf(fit)['pdf'] for fit in fits]

    #TODO: Try to emove the loop from here
    def produce_fits_name(self, fits):
        """NOTE: EXPERIMENTAL.
        Return a list with the ids of the fits"""
        return [fit.name for fit in fits]

    #TODO: Try to remove the loop from here
    def produce_fits_as(self, fits_pdf_config):
        """NOTE: EXPERIMENTAL. Return  the as value of the fits, reading
        it from the installed pdf"""
        return [pdf.AlphaS_MZ for pdf in fits_pdf_config]

    #TODO: Try to remove the loop from here.
    def produce_fits_as_from_fitdeclarations(self, fitdeclarations):
        """NOTE: EXPERIMENTAL. A hack to obtain fits_as from the
        fitdeclarations, without having to
        download and inspect the actual fits."""
        alpha_pattern = r'NNPDF\d\d_[a-z]+_as_(\d\d\d\d).*'
        res = []
        for fit in fitdeclarations:
            m = re.match(alpha_pattern, fit)
            if not m:
                raise ConfigError(f"Couldn't match fit name {fit} to the "
                                  "pattern {alpha_pattern!r}")
            res.append(float(m.group(1))/1000)
        return {'fits_as' : res}

    def produce_fits_name_from_fitdeclarations(self, fitdeclarations):
        """Inject the names from the ``fitdeclarations`` as the fit_names
        property"""
        #Cast nslist away
        return {'fits_name': list(fitdeclarations)}


    def produce_matched_datasets_from_dataspecs(self, dataspecs):
        """Take an arbitrary list of mappings called dataspecs and
        return a new list of mappings called dataspecs constructed as follows.

        From each of the original datasepcs, resolve the key `experiments` and
        all the dataset therein.

        Compute the intersection of the dataset names, and for each element in
        the intersection construct a mapping with the follwing keys:

            - experiment_name : A string with the common experiment name.
            - dataset_name : A string with the common dataset name.
            - datasepcs : A list of mappinngs matching the original
              "datasepcs". Each mapping contains:
                * dataset: A dataset with the name data_set name and the
                properties (cuts, theory, etc) corresponding to the original
                datasepec.
                * dataset_input: The input line used to build dataset.
                * All the other keys in the original dataspec.
        """
        if not isinstance(dataspecs, Sequence):
            raise ConfigError("dataspecs should be a sequence of mappings, not "
                              f"{type(dataspecs).__name__}")
        all_names = []
        for spec in dataspecs:
            if not isinstance(spec, Mapping):
                raise ConfigError("dataspecs should be a sequence of mappings, "
                      f" but {spec} is {type(spec).__name__}")

            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, experiments = self.parse_from_(None, 'experiments', write=False)
                names = {(e.name, ds.name):(ds, dsin) for e in experiments for ds, dsin in zip(e.datasets, e)}
                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))

        res = []
        for k in used_set:
            inres = {'experiment_name':k[0], 'dataset_name': k[1]}
            #TODO: Should this have the same name?
            l = inres['dataspecs'] = []
            for ispec, spec in enumerate(dataspecs):
                #Passing spec by referene
                d = ChainMap({
                    'dataset':       all_names[ispec][k][0],
                    'dataset_input': all_names[ispec][k][1],

                    },
                    spec)
                l.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x['experiment_name'], x['dataset_name']))
        return res

    #TODO: autogenerate functions like this
    def parse_experiments_covmat_output(self, fname:str, config_rel_path):
        """NOTE: THIS INTERFACE IS EXPERIMENTAL AND MIGHT CHANGE IN THE FUTURE.
        Process the output CSV table of the experiments_covmat action
        and return an equivalent datadrame"""
        from reportengine import filefinder

        finder = filefinder.FallbackFinder(['.', config_rel_path])
        try:
            folder, name = finder.find(fname)
        except FileNotFoundError as e:
            raise ConfigError(f"Couldn't locate '{fname}': {e}") from e
        try:
            df = tableloader.parse_experiments_covmat(folder/name)
        except Exception as e:
            raise ConfigError(e) from e

        return {'experiments_covmat': df}


    #TODO: Move these to their own module when that's supported by reportengine
    def produce_fits_matched_pseudorreplicas_chi2_output(self,
                                                         pseudorreplicafile:str,
                                                         fits):
        """NOTE: THIS IS A QUICK HACK AND MAY BE REMOVED OR HEAVILY CHANGED
        IN THE FUTURE (hence the horrible name).
        Quick and dirty analizing tool for the output of
        fits_matched_pseudorreplicas_chi2_table."""
        import pandas as pd
        import numpy as np
        try:
            df = pd.DataFrame.from_csv(pseudorreplicafile, sep='\t',
                index_col=[0,1],header=[0,1])
        except Exception as e:
            raise ConfigError(f"Failed to load the table: {e}") from e


        #Require that the fits are matched so we filer out some that are not
        #interesting or broken.
        try:
            df = df[[fit.name for fit in fits]]
        except Exception as e:
            raise ConfigError("Mismatch between fits provided and fits "
                             f"in the table {pseudorreplicafile}:\n{e}") from e
        ndataindexer = df.columns.get_locs([slice(None), 'ndata'])
        lentest = lambda x: len(np.unique(x.dropna()))<=1
        samelens = df.iloc[:,ndataindexer].apply(lentest, axis=1).all()
        if not samelens:
            raise ConfigError("Incorrect data: Expected all experiments to have the same length.")
        chindexer = df.columns.get_locs([slice(None), 'central_chi2'])
        df = df.iloc[:,chindexer]
        df = df.swaplevel(0,1)
        #Have it the way the existing functions like
        newcols = df.columns.set_levels([df.columns.levels[0], ['chi2']])
        df.columns = newcols
        return df

    def produce_fits_absolute_chi2_output(self, fitsexperimentsfile:str, fits):
        """"NOTE: THIS IS A QUICK HACK AND MAY BE REMOVED OR HEAVILY CHANGED
        IN THE FUTURE. Read the output of `fits_experiments_chi2_table` and
        produce a table with the total chi²."""
        import pandas as pd
        import numpy as np
        try:
            df = pd.DataFrame.from_csv(fitsexperimentsfile, sep='\t',
                header=[0,1], index_col=[0])
        except Exception as e:
            raise ConfigError(f"Failed to load the table: {e}") from e

        try:
            df = df[[fit.name for fit in fits]]
        except Exception as e:
            raise ConfigError("Mismatch between fits provided and fits "
                             f"in the table {fitsexperimentsfile}:\n{e}") from e

        #Multiply ndata * chi²/ndtata
        f = lambda x: x[x.columns[0]]*x[x.columns[1]]
        #Group by fit
        df = df.groupby(axis=1, level=0).apply(f)
        df.columns = pd.MultiIndex.from_product([np.array(df.columns), ['chi2']])
        return df

    def produce_fits_absolute_chi2_output_by_experiment(self,
            fits_absolute_chi2_output, prepend_total:bool=True):
        """Take the table returned by
        ``fits_absolute_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included."""
        if prepend_total:
            total = [{
            'fits_total_chi2': fits_absolute_chi2_output.sum(),
            'suptitle': 'Total'}]
        else:
            total = []

        return [*total,
                *[{'fits_total_chi2':df,
                 'suptitle': exp}
                for (exp, df) in
                fits_absolute_chi2_output.groupby(level=0)]
               ]




    def produce_fits_matched_pseudorreplicas_chi2_output_by_experiment(self,
            fits_matched_pseudorreplicas_chi2_output, prepend_total:bool=True):
        """Take the table returned by
        ``fits_matched_pseudorreplicas_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included.

        This provides a namespace list with `suptilte` and
        `fits_replica_data_correlated`.

        """

        if prepend_total:
            total = [{
            'fits_replica_data_correlated': fits_matched_pseudorreplicas_chi2_output.sum(level=1),
            'suptitle': 'Total'}]
        else:
            total = []


        return [*total,
                *[{'fits_replica_data_correlated': df,
                 'suptitle': exp}
                for (exp, df) in
                fits_matched_pseudorreplicas_chi2_output.groupby(level=0)]
               ]
