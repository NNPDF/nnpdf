%module(package="NNPDF") lhapdfset
 %{
#define SWIG_FILE_WITH_INIT
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/lhapdfset.h"
#include "LHAPDF/LHAPDF.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "include/numpy.i"

%init %{
    import_array();
%}

%import "pdfset.i"
%include "common.i"

%include "include/real_typemap.i"

%apply (NNPDF::real** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3)\
{(NNPDF::real** datamat, int* nrep_out, int* nf_out, int* nx_out)}

%apply (NNPDF::real* IN_ARRAY1, int DIM1) {(NNPDF::real* xmat, int nx_in)}
%apply (int* IN_ARRAY1, int DIM1) {(int* flmat, int nf_in)}


/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"
%include "include/real_typemap.i"

%include "../../src/NNPDF/lhapdfset.h"

namespace LHAPDF{
    int verbosity();
    void setVerbosity(int v);
}

%extend NNPDF::LHAPDFSet{

void grid_values(
                 int* flmat, int nf_in,
                 NNPDF::real* xmat, int nx_in,
                 NNPDF::real** datamat, int* nrep_out, int* nf_out, int* nx_out,
                 NNPDF::real Q
                )
{
    int nrep = self->GetMembers();
    *nrep_out = nrep;
    *nx_out = nx_in;
    *nf_out = nf_in;
    NNPDF::real* result = new NNPDF::real[nrep*nx_in*nf_in];
    for (int irep=0; irep<nrep; irep++){
        for (int ifl=0; ifl<nf_in; ifl++){
            for (int ix=0; ix<nx_in; ix++){
                int index = irep*nf_in*nx_in + ifl*nx_in + ix;
                result[index] =self->xfxQ(xmat[ix], Q, irep, flmat[ifl]);
            }
        }
    }
    *datamat = result;
}
}
