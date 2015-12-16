%module(package="NNPDF") nnmpi
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/nnmpi.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"
/* Parse the header file to generate wrappers */

%include "include/excepthandler.i"

%include "../../src/NNPDF/nnmpi.h"
