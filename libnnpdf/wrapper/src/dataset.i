%module(package="NNPDF") dataset
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/dataset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "commondata.i"
%import "fkset.i"
/* Parse the header file to generate wrappers */

%include "include/excepthandler.i"

%include "../../src/NNPDF/dataset.h"
