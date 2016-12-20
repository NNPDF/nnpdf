%module(package="NNPDF") nnmpi
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/nnmpi.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/nnmpi.h"
