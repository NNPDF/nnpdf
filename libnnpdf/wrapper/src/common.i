%module(package="NNPDF") common
 %{
#include "../../src/NNPDF/common.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/common.h"
