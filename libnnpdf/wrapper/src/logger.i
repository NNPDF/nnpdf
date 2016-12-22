%module(package="NNPDF") logger
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/logger.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/logger.h"
