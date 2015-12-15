%module(package="NNPDF") fkset
 %{
#include "../../src/NNPDF/fkset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%import "fastkernel.i"
%include "fkgenerator.i"
%include "../../src/NNPDF/fkset.h"
