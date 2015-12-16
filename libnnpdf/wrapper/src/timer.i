%module(package="NNPDF") timer
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/timer.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%include "include/excepthandler.i"

%include "../../src/NNPDF/timer.h"
