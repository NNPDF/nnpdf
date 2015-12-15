%module(package="NNPDF") fkgenerator
 %{
#include "../../src/NNPDF/fkgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "fastkernel.i"
%import "common.i"
/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/fkgenerator.h"
