%module(package="NNPDF") timer
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/timer.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/timer.h"
