%module(package="NNPDF") parametrisation
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/parametrisation.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 


/* Parse the header file to generate wrappers */
%import "common.i"
%include "include/excepthandler.i"

%include "../../src/NNPDF/parametrisation.h"
