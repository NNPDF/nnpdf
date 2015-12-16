%module(package="NNPDF") fkgenerator
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/fkgenerator.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "fastkernel.i"
%import "common.i"
/* Parse the header file to generate wrappers */

%include "include/excepthandler.i"

%include "../../src/NNPDF/fkgenerator.h"
