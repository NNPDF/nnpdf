%module(package="NNPDF") positivity
 %{
#define SWIG_FILE_WITH_INIT
#include "NNPDF/exceptions.h"
#include "NNPDF/positivity.h"
 %}
%include typemaps.i
%include "std_string.i" 
%include "std_vector.i" 


%include "include/numpy.i"

%init %{
    import_array();
%}

%include "commondata.i"
%include "fastkernel.i"
%include "common.i"


%include "include/real_typemap.i"

/* Parse the header file to generate wrappers */


%feature("autodoc", "3");

%include "include/excepthandler.i"

%apply int *OUTPUT {int *res};
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2)\
{(NNPDF::real** result, int* ndata, int* npdf)}

%include "NNPDF/positivity.h"
