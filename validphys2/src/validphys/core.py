# -*- coding: utf-8 -*-
"""
Core datastructures used in the validphys data model. Some of these are inmutable
specifications representing C++ objects.
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from __future__ import generator_stop

from collections import namedtuple
import functools
import inspect
import logging

import numpy as np
import scipy.stats

from reportengine import namespaces
from reportengine.baseexceptions import AsInputError

from NNPDF import LHAPDFSet
from NNPDF import CommonData, FKTable
from NNPDF.fkset import FKSet
from NNPDF.dataset import DataSet
from NNPDF.experiments import Experiment
from NNPDF.positivity import PositivitySet

from validphys import lhaindex

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
        argvals = ', '.join('%s=%r'%vals for vals in zip(self.argnames(),
                                                         self.comp_tuple))
        return '%s(%s)'%(self.__class__.__qualname__, argvals)

class PDFDoesNotExist(Exception): pass

class _PDFSETS():
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

    def __init__(self, name):
        self.name = name
        self._plotname = name
        super().__init__(name)


    def __getattr__(self, attr):
        #We don't even try to get reserved attributes from the info file
        if attr.startswith('__'):
            raise AttributeError(attr)
        try:
            return lhaindex.parse_info(self.name)[attr]
        except KeyError:
            raise AttributeError("'%r' has no attribute '%s'" % (type(self),
                                                                 attr))
        except IOError:
            raise PDFDoesNotExist(self.name)

    @property
    def label(self):
        return self._plotname

    @label.setter
    def label(self, label):
        self._plotname = label

    @property
    def stats_class(self):
        """Return the stats calculator for this error type"""
        error = self.ErrorType
        klass = STAT_TYPES[error]
        if hasattr(self, 'ErrorConfLevel'):
            klass = functools.partial(klass, rescale_factor=self.rescale_factor)
        return klass

    @property
    def infopath(self):
        return lhaindex.infofilename(self.name)

    @property
    def isinstalled(self):
        try:
            self.infopath
        except FileNotFoundError:
            return False
        else:
            return True


    @property
    def rescale_factor(self):
        if hasattr(self, "ErrorConfLevel"):
            if self.ErrorType == 'replicas':
                raise ValueError("Attribute at %s 'ErrorConfLevel' doesn't "
                "make sense with 'replicas' error type" % self.infopath)
            val = scipy.stats.norm.isf((1 - 0.01*self.ErrorConfLevel)/2)
            if np.isnan(val):
                raise ValueError("Invalid 'ErrorConfLevel' of PDF %s: %s" %
                                 (self, val))
            return val
        else:
            return 1

    @functools.lru_cache()
    def load(self):
        return LHAPDFSet(self.name, self.nnpdf_error)


    def __str__(self):
        return self.name

    def __len__(self):
        return self.NumMembers



    @property
    def nnpdf_error(self):
        """Return the NNPDF error tag, used to build the `LHAPDFSet` objeect"""
        error = self.ErrorType
        if error == "replicas":
            return LHAPDFSet.ER_MC
        if error == "hessian":
            if hasattr(self, 'ErrorConfLevel'):
                cl = self.ErrorConfLevel
                if cl == 90:
                    return LHAPDFSet.ER_EIG90
                elif cl == 68:
                    return LHAPDFSet.ER_EIG
                else:
                    raise NotImplementedError("No hessian errors with confidence"
                                              " interval %s" % (cl,) )
            else:
                return LHAPDFSet.ER_EIG

        if error == "symmhessian":
            if hasattr(self, 'ErrorConfLevel'):
                cl = self.ErrorConfLevel
                if cl == 68:
                    return LHAPDFSet.ER_SYMEIG
                else:
                    raise NotImplementedError("No symmetric hessian errors "
                                              "with confidence"
                                              " interval %s" % (cl,) )
            else:
                return LHAPDFSet.ER_EIG

        raise NotImplementedError("Error type for %s: '%s' is not implemented" %
                                  (self.name, error))

class CommonDataSpec(TupleComp):
    def __init__(self, datafile, sysfile, plotfiles):
        self.datafile = datafile
        self.sysfile = sysfile
        self.plotfiles = plotfiles
        super().__init__(datafile, sysfile, plotfiles)

    def __iter__(self):
        return iter((self.datafile, self.sysfile, self.plotfiles))

    @functools.lru_cache()
    def load(self):
        #TODO: Use better path handling in python 3.6
        return CommonData.ReadFile(str(self.datafile), str(self.sysfile))


class DataSetSpec(TupleComp):

    def __init__(self, *, name, commondata, fkspecs, thspec, cuts,
                 op=None):
        self.name = name
        self.commondata = commondata

        if isinstance(fkspecs, FKTableSpec):
            fkspecs = (fkspecs,)
        self.fkspecs = tuple(fkspecs)
        self.thspec = thspec

        self.cuts = cuts

        #Do this way (instead of setting op='NULL' in the signature)
        #so we don't have to know the default everywhere
        if op is None:
            op = 'NULL'
        self.op = op

        #TODO: Does it make sense to have a less verbose (but more obscure)
        #way to do this?
        #Note that we need to convert to tuple for hashing purposes
        if cuts is not None:
            frozencuts = cuts.data.tobytes()
        else:
            frozencuts = None
        super().__init__(name, commondata, fkspecs, thspec, frozencuts,
                         op)

    @functools.lru_cache()
    def load(self):
        cd = self.commondata.load()

        fktables = []
        for p in self.fkspecs:
            fktable = p.load()
            #IMPORTANT: We need to tell the python garbage collector to NOT free the
            #memory owned by the FKTable on garbage collection.
            #TODO: Do this automatically
            fktable.thisown = 0
            fktables.append(fktable)

        fkset = FKSet(FKSet.parseOperator(self.op), fktables)

        data = DataSet(cd, fkset)

        if self.cuts is not None:
            #ugly need to convert from numpy.int64 to int, so we can pass
            #it happily to the vector to the SWIG wrapper.
            #Do not do this (or find how to enable in SWIG):
            #data = DataSet(data, list(dataset.cuts))
            intmask = [int(ele) for ele in self.cuts]
            data = DataSet(data, intmask)
        return data

    def __str__(self):
        return self.name

class FKTableSpec(TupleComp):
    def __init__(self, fkpath, cfactors):
        self.fkpath = fkpath
        self.cfactors = cfactors
        super().__init__(fkpath, cfactors)

    #NOTE: We cannot do this because Fkset owns the fktable, and trying
    #to reuse the loaded one fails after it gets deleted.
    #@functools.lru_cache()
    def load(self):
        return FKTable(str(self.fkpath), [str(factor) for factor in self.cfactors])

class PositivitySetSpec(TupleComp):
    def __init__(self, commondataspec, fkspec, poslambda, thspec):
        self.commondataspec = commondataspec
        self.fkspec = fkspec
        self.poslambda = poslambda
        self.thspec = thspec
        super().__init__(commondataspec, fkspec, poslambda, thspec)

    @property
    def name(self):
        return self.commondataspec.name

    def __str__(self):
        return self.name

    @functools.lru_cache()
    def load(self):
        cd = self.commondataspec.load()
        fk = self.fkspec.load()
        return PositivitySet(cd, fk, self.poslambda)

    #__slots__ = ('__weakref__', 'commondataspec', 'fkspec', 'poslambda')


#We allow to expand the experiment as a list of datasets
class ExperimentSpec(TupleComp, namespaces.NSList):

    def __init__(self, name, datasets):
        #This needs to be hashable
        datasets = tuple(datasets)
        self.name = name
        self.datasets = datasets


        super().__init__(name, datasets)

        #TODO: Can we do  better cooperative inherece trick than this?
        namespaces.NSList.__init__(self, datasets, nskey='dataset')

    @functools.lru_cache()
    def load(self):
        sets = []
        for dataset in self.datasets:
            loaded_data = dataset.load()
            sets.append(loaded_data)
        return Experiment(sets, self.name)

    @property
    def thspec(self):
        #TODO: Is this good enough? Should we explicitly pass the theory
        return self.datasets[0].thspec


class FitSpec(TupleComp):
    def __init__(self, name, path):
        self.name = name
        self.path = path
        super().__init__(name, path)

    def __iter__(self):
        yield self.name
        yield self.path

    @functools.lru_cache()
    def as_input(self):
        import yaml
        p = self.path/'filter.yml'
        log.debug('Reading input from fit configuration %s' % (p,))
        try:
            with p.open() as f:
                d = yaml.load(f)
        except (yaml.YAMLError, FileNotFoundError) as e:
            raise AsInputError(str(e)) from e
        d['pdf'] = self.name
        return d

    __slots__ = ('label','name', 'path')



class TheoryIDSpec:
    def __init__(self, id, path):
        self.id = id
        self.path = path

    def __iter__(self):
        yield self.id
        yield self.path

    def __str__(self):
        return "theory %s" % self.id



#TODO: Decide if we want methods or properties
class Stats:

    def __init__(self, data):
        """`data `should be N_pdf*N_bins"""
        self.data = np.atleast_2d(data)

    def central_value(self):
        raise NotImplementedError()

    def error_members(self):
        raise NotImplementedError()

    def std_error(self):
        raise NotImplementedError()

    def sample_values(self, size):
        raise NotImplementedError()

    def std_interval(self, nsigma):
        raise NotImplementedError()

    def errorbar68(self):
        raise NotImplementedError()

    def errorbarstd(self):
        return (self.central_value() - self.std_error(),
                self.central_value() + self.std_error())

    #TODO...
    ...

class MCStats(Stats):
    """Result obtained from a Monte Carlo sample"""

    def central_value(self):
        return np.mean(self.data, axis=0)

    def std_error(self):
        return np.std(self.data, axis=0)

    def errorbar68(self):
        down = np.percentile(self.error_members(), 15.87, axis=0)
        up =   np.percentile(self.error_members(), 84.13, axis=0)
        return down, up

    def sample_values(self, size):
        return np.random.choice(self, size=size)

    def error_members(self):
        return self.data


class SymmHessianStats(Stats):
    """Compute stats in the 'assymetric' hessian format: The first index (0)
    is the
    central value. The rest of the indexes are results for each eigenvector.
    A 'rescale_factor is allowed in case the eigenvector confidence interval
    is not 68%'."""
    def __init__(self, data, rescale_factor=1):
        super().__init__(data)
        self.rescale_factor = rescale_factor

    def central_value(self):
        return self.data[0]

    def errorbar68(self):
        return self.errorbarstd()

    def error_members(self):
        return self.data[1:]

    def std_error(self):
        data = self.data
        diffsq = (data[0] - data[1:])**2
        return np.sqrt(diffsq.sum(axis=0))/self.rescale_factor

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
        diffsq = (data[1::2] - data[2::2])**2
        return np.sqrt(diffsq.sum(axis=0))/self.rescale_factor


STAT_TYPES = dict(
                    symmhessian = SymmHessianStats,
                    hessian = HessianStats,
                    replicas   = MCStats,
                   )

class Filter:
    def __init__(self, indexes, label, **kwargs):
        self.indexes  = indexes
        self.label = label
        self.kwargs = kwargs

    def as_pair(self):
        return self.label, self.indexes

    def __str__(self):
        return '%s: %s' % (self.label, self.indexes)
