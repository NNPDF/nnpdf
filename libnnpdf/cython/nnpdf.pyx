#cython embedsignature

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

from libcpp.string cimport string
from libcpp.vector cimport vector

#________________________________________________________
cdef extern from "lhapdfset.h":
    cdef cppclass _LHAPDFSet "NNPDF::LHAPDFSet":
        _LHAPDFSet(const string&, int) except +

cdef class LHAPDFSet:
    """\
    The PDF set from LHAPDF
    """
    cdef _LHAPDFSet *_ptr
    def __cinit__(self, const string& pdfname, int etype):
        self._ptr = new _LHAPDFSet(pdfname, etype)
    def __dealloc__(self):
        del self._ptr

class PDFSet:
    ER_NONE = 0
    ER_MC = 1
    ER_MC68 = 2
    ER_MCT0 = 3
    ER_EIG = 4
    ER_EIG90 = 5
    ER_SYMEIG = 6

#________________________________________________________
cdef extern from "commondata.h":
    cdef cppclass _CommonData "NNPDF::CommonData":
        _CommonData(const _CommonData& set)
        @staticmethod
        _CommonData ReadFile(const string&, const string&)

cdef class CommonData:
    """\
    The CommonData class
    """
    cdef _CommonData *_ptr
    cdef set_ptr(self, _CommonData* ptr):
        self._ptr = ptr

    @staticmethod
    def ReadFile(const string& filename, const string& systype):
        cdef _CommonData* ptr = new _CommonData(_CommonData.ReadFile(filename,systype))
        cdef CommonData obj
        obj = CommonData.__new__(CommonData)
        obj.set_ptr(ptr)
        return obj

#________________________________________________________
cdef extern from "fastkernel.h":
    cdef cppclass _FKTable "NNPDF::FKTable":
        _FKTable(const string&, const vector[string] &) except +
        int GetNData()

cdef class FKTable:
    """\
    The FKTable class
    """
    cdef _FKTable *_ptr
    def __cinit__(self, const string& filename, const vector[string]& cfactors):
        self._ptr = new _FKTable(filename, cfactors)
    def __dealloc__(self):
        del self._ptr
    @property
    def GetNData(self):
        return self._ptr.GetNData()

#________________________________________________________
cdef extern from "thpredictions.h":
    cdef cppclass _ThPredictions "NNPDF::ThPredictions":
        _ThPredictions(const _LHAPDFSet *, const _FKTable *) except +
        double GetObsCV(const int& idat)
        double GetObsError(const int& idat)
        int GetNData()
        string GetSetName()
        string GetPDFName()

cdef class ThPredictions:
    """\
    The ThPredictions class
    """
    cdef _ThPredictions *_ptr
    def __cinit__(self, LHAPDFSet pdfset, FKTable fktab):
        self._ptr = new _ThPredictions(pdfset._ptr, fktab._ptr)
    def __dealloc__(self):
        del self._ptr
    def GetObsCV(self, idat):
        return self._ptr.GetObsCV(idat)
    def GetObsError(self, idat):
        return self._ptr.GetObsError(idat)
    @property
    def GetNData(self):
        return self._ptr.GetNData()
    @property
    def GetSetName(self):
        return self._ptr.GetSetName()
    @property
    def GetPDFName(self):
        return self._ptr.GetPDFName()
