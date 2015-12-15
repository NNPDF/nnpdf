%module(package="NNPDF") dataset
 %{
#include "../../src/NNPDF/dataset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "commondata.i"
%import "fkset.i"
/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/dataset.h"
