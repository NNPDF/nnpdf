%module(package="NNPDF") randomgenerator
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/randomgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/randomgenerator.h"
