%module(package="NNPDF") logger
 %{
#include "../../src/NNPDF/logger.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/logger.h"
