%module(package="NNPDF") timer
 %{
#include "../../src/NNPDF/timer.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/timer.h"
