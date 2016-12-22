%module(package="NNPDF") parametrisation
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/parametrisation.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/parametrisation.h"
