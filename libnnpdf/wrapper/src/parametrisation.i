%module(package="NNPDF") parametrisation
 %{
#include "../../src/NNPDF/parametrisation.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */
%import "common.i"
%include "../../src/NNPDF/parametrisation.h"
