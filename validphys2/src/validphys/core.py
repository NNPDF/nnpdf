# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from collections import namedtuple
import functools

import numpy as np
import scipy.stats

from NNPDF import LHAPDFSet
from NNPDF import CommonData, FKTable
from NNPDF.fkset import FKSet
from NNPDF.dataset import DataSet


from validphys import lhaindex



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

class PDF:

    def __init__(self, name):
        self.name = name

    def __getattr__(self, attr):
        try:
            return lhaindex.parse_info(self.name)[attr]
        except KeyError:
            raise AttributeError("'%r' has no attribute '%s'" % (type(self),
                                                                 attr))
        except IOError:
            raise PDFDoesNotExist(self.name)

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

    def load(self):
        return LHAPDFSet(self.name, self.nnpdf_error)


    def __str__(self):
        return self.name



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


CommonDataSpec = namedtuple('CommonDataSpec', ['datafile', 'sysfile', 'plotfiles'])

class DataSetSpec:

    def __init__(self, *, name, commondata, cfac, fkpath, thspec, cuts):
        self.name = name
        self.commondata = commondata
        self.cfac = cfac
        self.fkpath = fkpath
        self.thspec = thspec

        self.cuts = cuts

    def load(self):
        cdpath,syspth, _ = self.commondata
        cd = CommonData.ReadFile(str(cdpath), str(syspth))

        fktable = FKTable(str(self.fkpath), [str(factor) for factor in self.cfac])
        #IMPORTANT: We need to tell the python garbage collector to NOT free the
        #memory owned by the FKTable on garbage collection.
        #TODO: Do this automatically
        fktable.thisown = 0
        fkset = FKSet(FKSet.parseOperator("NULL"), [fktable])

        data = DataSet(cd, fkset)

        if self.cuts is not None:
            #ugly need to convert from numpy.int64 to int, so we can pass
            #it happily to the vector to the SWIG wrapper.
            #Do not do this (or find how to enable in SWIG):
            #data = DataSet(data, list(dataset.cuts))
            intmask = [int(ele) for ele in self.cuts]
            data = DataSet(data, intmask)
        return data

FitSpec = namedtuple('FitSpec', ('name', 'path'))

class TheoryIDSpec:
    def __init__(self, id, path):
        self.id = id
        self.path = path

    def __iter__(self):
        yield self.id
        yield self.path


class Stats:

    def __init__(self, data):
        """`data `should be N_pdf*N_bins"""
        self.data = np.atleast_2d(data)

    def central_value(self):
        raise NotImplementedError()

    def std(self):
        raise NotImplementedError()

    def sample_values(self, size):
        raise NotImplementedError()

    def std_interval(self, nsigma):
        raise NotImplementedError()

    def errorbar68(self):
        raise NotImplementedError()

    #TODO...
    ...

class MCStats(Stats):
    """Result obtained from a Monte Carlo sample"""

    def central_value(self):
        return np.mean(self.data, axis=0)

    def std_error(self):
        return np.std(self.data, axis=0)

    def sample_values(self, size):
        return np.random.choice(self, size=size)


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
