%module(package="NNPDF") fkgenerator
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/fkgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "fastkernel.i"
%import "common.i"
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/fkgenerator.h"
