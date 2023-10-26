# -*- coding: utf-8 -*-
"""
Core datastructures used in the validphys data model.
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from __future__ import generator_stop

import enum
import functools
import inspect
import json
import logging
from pathlib import Path
import re

import numpy as np

from reportengine import namespaces
from reportengine.baseexceptions import AsInputError
from reportengine.compat import yaml

# TODO: There is a bit of a circular dependency between filters.py and this.
# Maybe move the cuts logic to its own module?
from validphys import filters, lhaindex
from validphys.commondataparser import (
    get_plot_kinlabels,
    parse_commondata,
    parse_commondata_new,
    peek_commondata_metadata,
)
from validphys.fkparser import load_fktable, parse_cfactor
from validphys.hyperoptplot import HyperoptTrial
from validphys.lhapdfset import LHAPDFSet
from validphys.tableloader import parse_exp_mat
from validphys.theorydbutils import fetch_theory
from validphys.utils import experiments_to_dataset_inputs

log = logging.getLogger(__name__)


class TupleComp:
    @classmethod
    def argnames(cls):
        return list(inspect.signature(cls.__init__).parameters.keys())[1:]

    def __init__(self, *args, **kwargs):
        self.comp_tuple = (*args, *kwargs.values())

    def __eq__(self, other):
        return self.comp_tuple == other.comp_tuple

    def __hash__(self):
        return hash(self.comp_tuple)

    def __repr__(self):
        argvals = ', '.join('%s=%r' % vals for vals in zip(self.argnames(), self.comp_tuple))
        return '%s(%s)' % (self.__class__.__qualname__, argvals)


class PDFDoesNotExist(Exception):
    pass


class _PDFSETS:
    """Convenient way to access installed PDFS, by e.g. tab completing
    in ipython."""

    def __getattr__(self, attr):
        if lhaindex.isinstalled(attr):
            return PDF(attr)
        raise AttributeError()

    def __dir__(self):
        return lhaindex.expand_local_names('*')


PDFSETS = _PDFSETS()


class PDF(TupleComp):
    """Base validphys PDF providing high level access to metadata.

    Statistical estimators which depends on the PDF type (MC, Hessian...)
    are exposed as a :py:class:`Stats` object through the :py:attr:`stats_class` attribute
    The LHAPDF metadata can directly be accessed through the :py:attr:`info` attribute


    Examples
    --------
    >>> from validphys.api import API
    >>> from validphys.convolution import predictions
    >>> args = {"dataset_input":{"dataset": "ATLASTTBARTOT"}, "theoryid":162, "use_cuts":"internal"}
    >>> ds = API.dataset(**args)
    >>> pdf = API.pdf(pdf="NNPDF40_nnlo_as_01180")
    >>> preds = predictions(ds, pdf)
    >>> preds.shape
    (3, 100)
    """

    def __init__(self, name):
        self.name = name
        self._plotname = name
        self._info = None
        self._stats_class = None
        super().__init__(name)

    @property
    def label(self):
        return self._plotname

    @label.setter
    def label(self, label):
        self._plotname = label

    @property
    def stats_class(self):
        """Return the stats calculator for this error type"""
        if self._stats_class is None:
            try:
                klass = STAT_TYPES[self.error_type]
            except KeyError:
                raise NotImplementedError(f"No Stats class for error type {self.error_type}")
            if self.error_conf_level is not None:
                klass = functools.partial(klass, rescale_factor=self._rescale_factor())
            self._stats_class = klass
        return self._stats_class

    @property
    def infopath(self):
        return Path(lhaindex.infofilename(self.name))

    @property
    def info(self):
        """Information contained in the LHAPDF .info file"""
        if self._info is None:
            try:
                self._info = lhaindex.parse_info(self.name)
            except IOError as e:
                raise PDFDoesNotExist(self.name) from e
        return self._info

    @property
    def q_min(self):
        """Minimum Q as given by the LHAPDF .info file"""
        return self.info["QMin"]

    @property
    def error_type(self):
        """Error type as defined in the LHAPDF .info file"""
        return self.info["ErrorType"]

    @property
    def alphas_mz(self):
        """Alpha_s(M_Z) as defined in the LHAPDF .info file"""
        return self.info["AlphaS_MZ"]

    @property
    def alphas_vals(self):
        """List of alpha_s(Q) at various Q for interpolation based alphas.
        Values as defined in the LHAPDF .info file"""
        self.info["AlphaS_Vals"]

    @property
    def error_conf_level(self):
        """Error confidence level as defined in the LHAPDF .info file
        if no number is given in the LHAPDF .info file defaults to 68%
        """
        key_name = "ErrorConfLevel"
        # Check possible misconfigured info file
        if self.error_type == "replicas":
            if key_name in self.info:
                raise ValueError(
                    f"Attribute at {self.infopath} 'ErrorConfLevel' doesn't "
                    "make sense with 'replicas' error type"
                )
            return None
        return self.info.get(key_name, 68)

    @property
    def isinstalled(self):
        try:
            return self.infopath.is_file()
        except FileNotFoundError:
            return False

    def _rescale_factor(self):
        """Compute the rescale factor for the stats class"""
        # This is imported here for performance reasons.
        import scipy.stats

        val = scipy.stats.norm.isf((1 - 0.01 * self.error_conf_level) / 2)
        if np.isnan(val):
            raise ValueError(f"Invalid 'ErrorConfLevel' for PDF {self}: {val}")
        return val

    @functools.lru_cache(maxsize=16)
    def load(self):
        return LHAPDFSet(self.name, self.error_type)

    @functools.lru_cache(maxsize=2)
    def load_t0(self):
        """Load the PDF as a t0 set"""
        return LHAPDFSet(self.name, "t0")

    def __str__(self):
        return self.label

    def __len__(self):
        return self.info["NumMembers"]

    def get_members(self):
        """Return the number of members selected in ``pdf.load().grid_values``"""
        return len(self)


class CommonDataSpec(TupleComp):
    def __init__(self, datafile, sysfile, plotfiles, name=None, metadata=None, legacy=True):
        self.legacy = legacy
        self._metadata = metadata
        self.datafile = datafile
        self.sysfile = sysfile
        if legacy:
            self.plotfiles = tuple(plotfiles)
            super().__init__(datafile, sysfile, self.plotfiles)
        else:
            self.plotfiles = False
            super().__init__(name)

    @property
    def name(self):
        return self.metadata.name

    @property
    def nsys(self):
        if self.legacy:
            return self.metadata.nsys
        else:
            cd = self.load()
            return cd.nsys

    @property
    def ndata(self):
        if self.legacy:
            return self.metadata.ndata
        else:
            cd = self.load()
            return cd.nsys

    @property
    def process_type(self):
        return self.metadata.process_type

    @property
    def metadata(self):
        if self._metadata is None:
            self._metadata = peek_commondata_metadata(self.datafile)
        return self._metadata

    def __str__(self):
        return self.name

    def __iter__(self):
        return iter((self.datafile, self.sysfile, self.plotfiles))

    # TODO: one of the two functions below needs to go
    @functools.lru_cache()
    def load(self):
        if self.legacy:
            return parse_commondata(self.datafile, self.sysfile, self.name)
        else:
            return parse_commondata_new(self.metadata)

    def load_commondata_instance(self):
        """
        load a validphys.core.CommonDataSpec to validphys.core.CommonData
        """
        from validphys.commondataparser import load_commondata

        return load_commondata(self)

    @property
    def plot_kinlabels(self):
        if self.legacy:
            return get_plot_kinlabels(self)
        else:
            return self.metadata.kinlabels


class DataSetInput(TupleComp):
    """Represents whatever the user enters in the YAML to specify a
    dataset."""

    def __init__(self, *, name, sys, cfac, frac, weight, custom_group, variants):
        self.name = name
        self.sys = sys
        self.cfac = cfac
        self.frac = frac
        self.weight = weight
        self.custom_group = custom_group
        self.variants = variants
        super().__init__(name, sys, cfac, frac, weight, custom_group)

    def __str__(self):
        return self.name


class ExperimentInput(TupleComp):
    def __init__(self, *, name, datasets):
        self.name = name
        self.datasets = datasets
        super().__init__(name, datasets)

    def as_dict(self):
        return {'experiment': self.name, 'datasets': self.datasets}

    def __str__(self):
        return self.name


class CutsPolicy(enum.Enum):
    INTERNAL = "internal"
    NOCUTS = "nocuts"
    FROMFIT = "fromfit"
    FROM_CUT_INTERSECTION_NAMESPACE = "fromintersection"
    FROM_SIMILAR_PREDICTIONS_NAMESPACE = "fromsimilarpredictions"


class Cuts(TupleComp):
    def __init__(self, commondata, path):
        """Represents a file containing cuts for a given dataset"""
        name = commondata.name
        self.name = name
        self.path = path
        self._mask = None
        if path is None:
            log.debug("No filter found for %s, all points allowed", name)
            self._mask = np.arange(commondata.ndata)
        self._legacy_mask = None
        super().__init__(name, path)

    def load(self):
        if self._mask is not None:
            return self._mask
        log.debug("Loading cuts for %s", self.name)
        return np.atleast_1d(np.loadtxt(self.path, dtype=int))


class InternalCutsWrapper(TupleComp):
    def __init__(self, commondata, rules):
        self.rules = rules
        self.commondata = commondata
        super().__init__(commondata, tuple(rules))

    def load(self):
        return np.atleast_1d(
            np.asarray(filters.get_cuts_for_dataset(self.commondata, self.rules), dtype=int)
        )


class MatchedCuts(TupleComp):
    def __init__(self, othercuts, ndata):
        self.othercuts = tuple(othercuts)
        self.ndata = ndata
        super().__init__(self.othercuts, self.ndata)

    def load(self):
        loaded = [c.load() for c in self.othercuts if c]
        if loaded:
            return functools.reduce(np.intersect1d, loaded)
        self._full = True
        return np.arange(self.ndata)


class SimilarCuts(TupleComp):
    def __init__(self, inputs, threshold):
        if len(inputs) != 2:
            raise ValueError("Expecting two input tuples")
        firstcuts, secondcuts = inputs[0][0].cuts, inputs[1][0].cuts
        if firstcuts != secondcuts:
            raise ValueError("Expecting cuts to be the same for all datasets")
        self.inputs = inputs
        self.threshold = threshold
        super().__init__(self.inputs, self.threshold)

    @functools.lru_cache()
    def load(self):
        # TODO: Update this when a suitable interace becomes available
        from validphys.commondataparser import load_commondata
        from validphys.convolution import central_predictions
        from validphys.covmats import covmat_from_systematics

        first, second = self.inputs
        first_ds = first[0]
        exp_err = np.sqrt(
            np.diag(
                covmat_from_systematics(
                    load_commondata(first_ds.commondata).with_cuts(first_ds.cuts),
                    first_ds,  # DataSetSpec has weight attr
                    use_weights_in_covmat=False,  # Don't weight covmat
                )
            )
        )
        # Compute matched predictions
        delta = np.abs((central_predictions(*first) - central_predictions(*second)).squeeze(axis=1))
        ratio = delta / exp_err
        passed = ratio < self.threshold
        return passed[passed].index


def cut_mask(cuts):
    """Return an objects that will act as the cuts when applied as a slice"""
    if cuts is None:
        return slice(None)
    return cuts.load()


class DataSetSpec(TupleComp):
    def __init__(self, *, name, commondata, fkspecs, thspec, cuts, frac=1, op=None, weight=1):
        self.name = name
        self.commondata = commondata

        if isinstance(fkspecs, FKTableSpec):
            fkspecs = (fkspecs,)
        fkspecs = tuple(fkspecs)
        self.fkspecs = fkspecs
        self.thspec = thspec

        self.cuts = cuts
        self.frac = frac

        # Do this way (instead of setting op='NULL' in the signature)
        # so we don't have to know the default everywhere
        if op is None:
            op = 'NULL'
        self.op = op
        self.weight = weight

        super().__init__(name, commondata, fkspecs, thspec, cuts, frac, op, weight)

    @functools.lru_cache()
    def load_commondata(self):
        """Strips the commondata loading from `load`"""
        cd = self.commondata.load()
        if self.cuts is not None:
            loaded_cuts = self.cuts.load()
            if not (hasattr(loaded_cuts, '_full') and loaded_cuts._full):
                intmask = [int(ele) for ele in loaded_cuts]
                cd = cd.with_cuts(intmask)
        return cd

    def to_unweighted(self):
        """Return a copy of the dataset with the weight set to one."""
        return self.__class__(
            name=self.name,
            commondata=self.commondata,
            fkspecs=self.fkspecs,
            thspec=self.thspec,
            cuts=self.cuts,
            frac=self.frac,
            op=self.op,
            weight=1,
        )

    def __str__(self):
        return self.name


class FKTableSpec(TupleComp):
    """
    Each FKTable is formed by a number of sub-fktables to be concatenated
    each of which having its own path.
    Therefore the ``fkpath`` variable is a list of paths.

    Before the pineappl implementation, FKTable were already pre-concatenated.
    The Legacy interface therefore relies on fkpath being just a string or path instead

    The metadata of the FKTable for the given dataset is stored as an attribute to this function.
    This is transitional, eventually it will be held by the associated CommonData in the new format.
    """

    def __init__(self, fkpath, cfactors, metadata=None):
        self.cfactors = cfactors if cfactors is not None else []

        self.legacy = False

        # NOTE: The legacy interface is currently used by fkparser to decide
        # whether to read an FKTable using the old parser or the pineappl parser
        # this attribute (and the difference between both) might be removed in future
        # releases of NNPDF so please don't write code that relies on it
        if not isinstance(fkpath, (tuple, list)):
            self.legacy = True
        else:
            # Make it into a tuple only for the new format
            fkpath = tuple(fkpath)

        self.fkpath = fkpath
        self.metadata = metadata

        # For new theories, add also the target_dataset so that we don't reuse fktables
        # Ideally this won't be necessary in the future and we will be able to reutilize fktables.
        if not self.legacy:
            super().__init__(fkpath, cfactors, self.metadata.target_dataset)
        else:
            super().__init__(fkpath, cfactors)

    def load_cfactors(self):
        """Each of the sub-fktables that form the complete FKTable can have several cfactors
        applied to it. This function uses ``parse_cfactor`` to make them into CFactorData
        """
        if self.legacy:
            raise NotImplementedError("cfactor loading from spec not implemented for old theories")

        return [[parse_cfactor(c.open("rb")) for c in cfacs] for cfacs in self.cfactors]

    def load_with_cuts(self, cuts):
        """Load the fktable and apply cuts inmediately. Returns a FKTableData"""
        return load_fktable(self).with_cuts(cuts)


class LagrangeSetSpec(DataSetSpec):
    """Extends DataSetSpec to work around the particularities of the positivity, integrability
    and other Lagrange Multiplier datasets.
    """

    def __init__(self, name, commondataspec, fkspec, maxlambda, thspec):
        cuts = Cuts(commondataspec, None)
        self.maxlambda = maxlambda
        super().__init__(
            name=name,
            commondata=commondataspec,
            fkspecs=fkspec,
            thspec=thspec,
            cuts=cuts,
        )

    def to_unweighted(self):
        log.warning(
            "Trying to unweight %s, %s are always unweighted",
            self.__class__.__name__,
            self.name,
        )
        return self

    @functools.lru_cache()
    def load_commondata(self):
        return self.commondata.load()


class PositivitySetSpec(LagrangeSetSpec):
    pass


class IntegrabilitySetSpec(LagrangeSetSpec):
    pass


# We allow to expand the experiment as a list of datasets
class DataGroupSpec(TupleComp, namespaces.NSList):
    def __init__(self, name, datasets, dsinputs=None):
        # This needs to be hashable
        datasets = tuple(datasets)

        # TODO: Find a better way for interactive usage.
        if dsinputs is not None:
            dsinputs = tuple(dsinputs)

        self.name = name
        self.datasets = datasets
        self.dsinputs = dsinputs

        # TODO: Add dsinputs to comp tuple?
        super().__init__(name, datasets)

        # TODO: Can we do  better cooperative inherece trick than this?
        namespaces.NSList.__init__(self, dsinputs, nskey='dataset_input')

    @functools.lru_cache(maxsize=32)
    def load_commondata(self):
        return [d.load_commondata() for d in self.datasets]

    def load_commondata_instance(self):
        """
        Given Experiment load list of validphys.coredata.CommonData
        objects with cuts already applied
        """
        commodata_list = []
        for dataset in self.datasets:
            cd = dataset.commondata.load_commondata_instance()
            if dataset.cuts is None:
                commodata_list.append(cd)
            else:
                commodata_list.append(cd.with_cuts(dataset.cuts.load()))
        return commodata_list

    @property
    def thspec(self):
        # TODO: Is this good enough? Should we explicitly pass the theory
        return self.datasets[0].thspec

    def __str__(self):
        return self.name

    # Need this so that it doesn't try to iterte over itself.
    @property
    def as_markdown(self):
        return str(self)

    def to_unweighted(self):
        """Return a copy of the group with the weights for all experiments set
        to one. Note that the results cannot be used as a namespace."""
        return self.__class__(
            name=self.name,
            datasets=[ds.to_unweighted() for ds in self.datasets],
            dsinputs=None,
        )


class FitSpec(TupleComp):
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.label = name
        super().__init__(name, path)

    def __iter__(self):
        yield self.name
        yield self.path

    @functools.lru_cache()
    def as_input(self):
        p = self.path / 'filter.yml'
        log.debug('Reading input from fit configuration %s', p)
        try:
            with p.open() as f:
                d = yaml.safe_load(f)
        except (yaml.YAMLError, FileNotFoundError) as e:
            raise AsInputError(str(e)) from e
        d['pdf'] = {'id': self.name, 'label': self.label}

        if 'experiments' in d:
            # Flatten old style experiments to dataset_inputs
            dataset_inputs = experiments_to_dataset_inputs(d['experiments'])
            d['dataset_inputs'] = dataset_inputs

        # BCH
        # backwards compatibility hack for runcards with the 'fitting' namespace
        # if a variable already exists outside 'fitting' it takes precedence
        fitting = d.get("fitting")
        if fitting is not None:
            to_take_out = ["genrep", "trvlseed", "mcseed", "nnseed", "parameters"]
            for vari in to_take_out:
                if vari in fitting and vari not in d:
                    d[vari] = fitting[vari]
        return d

    def __str__(self):
        return self.label

    __slots__ = ('label', 'name', 'path')


class HyperscanSpec(FitSpec):
    """The hyperscan spec is just a special case of FitSpec"""

    def __init__(self, name, path):
        super().__init__(name, path)
        self._tries_files = None

    @property
    def tries_files(self):
        """Return a dictionary with all tries.json files mapped to their replica number"""
        if self._tries_files is None:
            re_idx = re.compile(r"(?<=replica_)\d+$")
            get_idx = lambda x: int(re_idx.findall(x.as_posix())[-1])
            all_rep = map(get_idx, self.path.glob("nnfit/replica_*"))
            # Now loop over all replicas and save them when they include a tries.json file
            tries = {}
            for idx in sorted(all_rep):
                test_path = self.path / f"nnfit/replica_{idx}/tries.json"
                if test_path.exists():
                    tries[idx] = test_path
            self._tries_files = tries
        return self._tries_files

    def get_all_trials(self, base_params=None):
        """Read all trials from all tries files.
        If there are original runcard-based parameters, a reference to them can be passed
        to the trials so that a full hyperparameter dictionary can be defined

        Each hyperopt trial object will also have a reference to all trials in its own file
        """
        all_trials = []
        for trial_file in self.tries_files.values():
            with open(trial_file, "r") as tf:
                run_trials = []
                for trial in json.load(tf):
                    trial = HyperoptTrial(trial, base_params=base_params, linked_trials=run_trials)
                    run_trials.append(trial)
            all_trials += run_trials
        return all_trials

    def sample_trials(self, n=None, base_params=None, sigma=4.0):
        """Parse all trials in the hyperscan object
        and then return an array of ``n`` trials read from the ``tries.json`` files
        and sampled according to their reward.
        If ``n`` is ``None``, no sapling is performed and all trials are returned

        Returns
        -------
            Dictionary on the form {parameters: list of trials}
        """
        all_trials_raw = self.get_all_trials(base_params=base_params)
        # Drop all failing trials
        all_trials = list(filter(lambda i: i.reward, all_trials_raw))
        if n is None:
            return all_trials
        if len(all_trials) < n:
            log.warning("Asked for %d trials, only %d valid trials found", n, len(all_trials))
        # Compute weights proportionally to the reward (goes from 0 (worst) to 1 (best, loss=1))
        rewards = np.array([i.weighted_reward for i in all_trials])
        weight_raw = np.exp(sigma * rewards**2)
        total = np.sum(weight_raw)
        weights = weight_raw / total
        return np.random.choice(all_trials, replace=False, size=n, p=weights)


class TheoryIDSpec:
    def __init__(self, id, path):
        self.id = id
        self.path = path

    def __iter__(self):
        yield self.id
        yield self.path

    def get_description(self):
        dbpath = self.path.parent / 'theory.db'
        return fetch_theory(dbpath, self.id)

    __slots__ = ('id', 'path')

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id}, path={self.path!r})"

    def __str__(self):
        return f"Theory {self.id}"

    @property
    def yamldb_path(self):
        return self.path / "yamldb"

    def is_pineappl(self):
        """Check whether this theory is a pineappl-based theory"""
        return self.yamldb_path.exists()


class ThCovMatSpec:
    def __init__(self, path):
        self.path = path

    # maxsize relatively low here, expect single experiments so one load per dataspec
    @functools.lru_cache(maxsize=8)
    def load(self):
        return parse_exp_mat(self.path)

    def __str__(self):
        return str(self.path)


# TODO: Decide if we want methods or properties
class Stats:
    def __init__(self, data):
        """`data `should be N_pdf*N_bins"""
        self.data = np.atleast_2d(data)

    def central_value(self):
        return self.data[0]

    def error_members(self):
        return self.data[1:]

    def std_error(self):
        raise NotImplementedError()

    def moment(self, order):
        raise NotImplementedError()

    def sample_values(self, size):
        raise NotImplementedError()

    def std_interval(self, nsigma):
        raise NotImplementedError()

    def errorbar68(self):
        raise NotImplementedError()

    def errorbarstd(self):
        return (self.central_value() - self.std_error(), self.central_value() + self.std_error())

    # TODO...
    ...


class MCStats(Stats):
    """Result obtained from a Monte Carlo sample"""

    def std_error(self):
        # ddof == 1 to match legacy libNNPDF behaviour
        return np.std(self.error_members(), ddof=1, axis=0)

    def moment(self, order):
        return np.mean(np.power(self.error_members() - self.central_value(), order), axis=0)

    def errorbar68(self):
        # Use nanpercentile here because we can have e.g. 0/0==nan normalization
        # somewhere.
        down = np.nanpercentile(self.error_members(), 15.87, axis=0)
        up = np.nanpercentile(self.error_members(), 84.13, axis=0)
        return down, up

    def sample_values(self, size):
        return np.random.choice(self, size=size)


class SymmHessianStats(Stats):
    """Compute stats in the 'symetric' hessian format: The first index (0)
    is the
    central value. The rest of the indexes are results for each eigenvector.
    A 'rescale_factor is allowed in case the eigenvector confidence interval
    is not 68%'."""

    def __init__(self, data, rescale_factor=1):
        super().__init__(data)
        self.rescale_factor = rescale_factor

    def errorbar68(self):
        return self.errorbarstd()

    def std_error(self):
        data = self.data
        diffsq = (data[0] - data[1:]) ** 2
        return np.sqrt(diffsq.sum(axis=0)) / self.rescale_factor

    def moment(self, order):
        data = self.data
        return np.sum(np.power((data[0] - data[1:]) / self.rescale_factor, order), axis=0)


class HessianStats(SymmHessianStats):
    """Compute stats in the 'assymetric' hessian format: The first index (0)
    is the
    central value. The odd indexes are the
    results for lower eigenvectors and the
    even are the upper eigenvectors.A 'rescale_factor is allowed in
    case the eigenvector confidence interval
    is not 68%'."""

    def std_error(self):
        data = self.data
        diffsq = (data[1::2] - data[2::2]) ** 2
        return np.sqrt(diffsq.sum(axis=0)) / self.rescale_factor / 2

    def moment(self, order):
        data = self.data
        return np.sum(np.power((data[1::2] - data[2::2]) / self.rescale_factor / 2, order), axis=0)


STAT_TYPES = dict(
    symmhessian=SymmHessianStats,
    hessian=HessianStats,
    replicas=MCStats,
)


class Filter:
    def __init__(self, indexes, label, **kwargs):
        self.indexes = indexes
        self.label = label
        self.kwargs = kwargs

    def as_pair(self):
        return self.label, self.indexes

    def __str__(self):
        return '%s: %s' % (self.label, self.indexes)
