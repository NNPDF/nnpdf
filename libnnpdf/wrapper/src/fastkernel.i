%module(package="NNPDF") fastkernel
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/fastkernel.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */
%include "NNPDF/exceptions.h"

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/fastkernel.h"
