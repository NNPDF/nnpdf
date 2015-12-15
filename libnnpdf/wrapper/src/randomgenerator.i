%module(package="NNPDF") randomgenerator
 %{
#include "../../src/NNPDF/randomgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/randomgenerator.h"
