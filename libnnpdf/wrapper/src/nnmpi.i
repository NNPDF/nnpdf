%module(package="NNPDF") nnmpi
 %{
#include "../../src/NNPDF/nnmpi.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"
/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/nnmpi.h"
