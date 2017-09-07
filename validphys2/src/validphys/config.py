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
        except LoaderError as e:
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
        known_keys = {'dataset', 'sys', 'cfac', 'frac'}
        try:
            name = dataset['dataset']
            if not isinstance(name, str):
                raise Config("'name' must be a string, not %s" % type(name))
        except KeyError:
            raise ConfigError("'dataset' must be a mapping with "
                              "'dataset' and 'sysnum'")


        sysnum = dataset.get('sys')
        cfac = dataset.get('cfac', tuple())
        kdiff = dataset.keys() - known_keys
        for k in kdiff:
            #Abuse ConfigError to get the suggestions.
            log.warn(ConfigError(f"Key '{k}' in dataset_input not known.", k, known_keys))
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

    def parse_blacklist_datasets(self, datasets:list):
        return datasets

    def produce_combine_dataspecs_pseudorreplicas_as(
            self, dataspecs, how='min', blacklist_datasets=()):
        if not isinstance(dataspecs, Sequence):
            raise ConfigError("dataspecs should be a sequence of mappings, not "
                              f"{type(dataspecs).__name__}")
        if how != 'min':
            raise ConfigError("Only min is implemented at the moment")

        dfs = []
        fitnames = []
        for spec in dataspecs:
            if not isinstance(spec, Mapping):
                raise ConfigError("dataspecs should be a sequence of mappings, "
                      f" but {spec} is {type(spec).__name__}")
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, df = self.parse_from_(None, 'fits_computed_psedorreplicas_chi2', write=False)
                _, asval = self.parse_from_(None, 'fits_as', write=False)
                _, namelist = self.parse_from_(None, 'fits_name', write=False)
                if not dfs:
                    firstas = asval
                elif asval != firstas:
                    raise ConfigError("Expecting all as values to be the same")
                dfs.append(df)
                fitnames.append(namelist)
        finalnames =  [min(ns, key=len) + '__combined' for ns in zip(*fitnames)]
        res = tableloader.combine_pseudorreplica_tables(dfs, finalnames,
                blacklist_datasets=blacklist_datasets)

        return {'fits_computed_psedorreplicas_chi2': res}



    def _get_table(self, loader_func, fname, config_rel_path):
        try:
            res = self.loader.check_vp_output_file(fname.strip(),
                    extra_paths=['.', config_rel_path])
        except LoaderError as e:
            raise ConfigError(e) from e

        try:
            df = loader_func(res)
        except Exception as e:
            raise ConfigError(e) from e
        return df


    #TODO: autogenerate functions like this
    def parse_experiments_covmat_output(self, fname:str, config_rel_path):
        """NOTE: THIS INTERFACE IS EXPERIMENTAL AND MIGHT CHANGE IN THE FUTURE.
        Process the output CSV table of the experiments_covmat action
        and return an equivalent datadrame"""
        df = self._get_table(tableloader.load_experiments_covmat, fname, config_rel_path)
        return {'experiments_covmat': df}


    #TODO: Move these to their own module when that's supported by reportengine
    def produce_fits_matched_pseudorreplicas_chi2_output(self,
                                                         pseudorreplicafile:str,
                                                         fits_name):
        """DEPRECATED. DO NOT USE."""
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
            df = df[fits_name]
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

    def parse_fits_computed_psedorreplicas_chi2_output(self, fname:str,
            config_rel_path):
        """Return a namespace (mapping) with the output of
        ``fits_computed_psedorreplicas_chi2_table`` as read from the specified
        filename. Use a {@with@} block to pass it to the providers.
        The fit names must be provided explicitly."""
        return self._get_table(tableloader.load_fits_computed_psedorreplicas_chi2,
                             fname, config_rel_path)


    def produce_use_fits_computed_psedorreplicas_chi2_output(
            self, fits_computed_psedorreplicas_chi2_output, fits_name):
        """Select the columns of the input file matching the fits."""
        df = fits_computed_psedorreplicas_chi2_output
        try:
            df = df[fits_name]
        except Exception as e:
            raise ConfigError(f"Could not select the fit names from the table: {e}") from e

        return {'fits_computed_psedorreplicas_chi2':  df}


    @element_of('extra_sums')
    def parse_extra_sum(self, s:dict):
        keys = {'dataset_item', 'components'}
        if s.keys() != keys:
            d1 = s.keys() - keys
            d2 = keys - s.keys
            if d1:
                raise ConfigError(f'Unable to parse extra_sum: unrecognized keys: {d1}')
            if d2:
                raise ConfigError(f'Unable to parse extra_sum. The following keys are required: {d2}')
            raise RuntimeError()
        return s


    def produce_fits_matched_pseudorreplicas_chi2_by_experiment_and_dataset(self,
            fits_computed_psedorreplicas_chi2, prepend_total:bool=True,
            extra_sums=None):
        """Take the table returned by
        ``fits_matched_pseudorreplicas_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included.

        This provides a namespace list with `suptilte`, `ndata` and
        `fits_replica_data_correlated`.

        """
        df = fits_computed_psedorreplicas_chi2

        if prepend_total:
            s =  df.loc[(slice(None), 'Total'),:].groupby(level=3).sum()
            ndata = sum(n for (n,_) in df.loc[(slice(None), 'Total'),:].groupby(level=2))
            total = [
                {'experiment_label': 'Total',
                'by_dataset': [{
                    'fits_replica_data_correlated': s,
                    'suptitle': 'Total',
                    'ndata': ndata
                 }]}]
        else:
            total = []

        expres = []
        for exp, expdf in df.groupby(level=0):
            d = {'experiment_label': exp}
            by_dataset = d['by_dataset'] = []
            for ds, dsdf in expdf.groupby(level=1):
                ndata = sum(n for (n,_) in dsdf.groupby(level=2))
                dsdf.index  = dsdf.index.droplevel([0,1,2])

                if ds == 'Total':
                    if exp != 'Total':
                        ds = f'{exp} Total'
                    by_dataset.insert(0, {'fits_replica_data_correlated': dsdf,
                                   'suptitle':ds, 'ndata':ndata})
                else:
                    by_dataset.append({'fits_replica_data_correlated': dsdf,
                                   'suptitle':ds, 'ndata':ndata})

            expres.append(d)

        if extra_sums:
            for es in extra_sums:
                label = es['dataset_item']
                components = es['components']
                dss = set(df.index.levels[1])
                diff = set(components) - dss
                if diff:
                    bad_item = next(iter(diff))
                    raise ConfigError(f"Unrecognized elements in extra_sum: {diff}", bad_item, dss)

                sliced = df.loc[(slice(None), components),:]
                s =  sliced.groupby(level=3).sum()
                ndata = sum(n for (n,_) in sliced.groupby(level=2))
                total.append(
                    {'experiment_label': label,
                    'by_dataset': [{
                        'fits_replica_data_correlated': s,
                        'suptitle': label,
                        'ndata': ndata
                     }]})


        return [*total, *expres]

    def parse_fits_chi2_paramfits_output(self, fname:str, config_rel_path):
        """Load the output of ``fits_chi2_table`` adapted to suit the
        ``paramfits`` module. The fit names must be provided explicitly."""
        return self._get_table(tableloader.load_adapted_fits_chi2_table,
                             fname, config_rel_path)

    def produce_use_fits_chi2_paramfits_output(self, fits_chi2_paramfits_output,
                                                   fits_name):
        df = fits_chi2_paramfits_output
        try:
            df = df[fits_name]
        except Exception as e:
            raise ConfigError(f"Could not select the fit names from the table: {e}") from e
        return {'adapted_fits_chi2_table':  df}

    def produce_fits_central_chi2_by_experiment_and_dataset(self,
            adapted_fits_chi2_table, prepend_total=True):
        """Take the table returned by
        ``fits_matched_pseudorreplicas_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included.

        This provides a namespace list with `suptilte` and
        `fits_replica_data_correlated`."""

        df = adapted_fits_chi2_table

        if prepend_total:
            s =  df.sort_index().loc[(slice(None), 'Total'), :].sum()
            total = [
                {'experiment_label': 'Total',
                'by_dataset': [{
                    'fits_total_chi2': s,
                    'suptitle': 'Total',
                 }]}]
        else:
            total = []
        expres = []
        for exp, expdf in df.groupby(level='experiment'):
            d = {'experiment_label': exp}
            by_dataset = d['by_dataset'] = []
            for ds, dsdf in expdf.groupby(level=1):
                dsdf.index  = dsdf.index.droplevel([0])
                if ds == 'Total':
                    ds = f'{exp} Total'
                    by_dataset.insert(0, {'fits_total_chi2': dsdf,
                                   'suptitle':ds})
                else:
                    by_dataset.append({'fits_total_chi2': dsdf,
                                   'suptitle':ds})

            expres.append(d)

        return [*total, *expres]



