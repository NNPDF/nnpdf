%module(package="NNPDF") randomgenerator
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/randomgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/randomgenerator.h"
