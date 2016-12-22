%module(package="NNPDF") common
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/common.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/common.h"

/* We define often used templates here to avoid name conflicts */
%template(vector_str) std::vector<std::string>;

