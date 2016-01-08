%module(package="NNPDF") logger
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/logger.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/logger.h"
