%module(package="NNPDF") fastkernel
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/fastkernel.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/exceptions.h"

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/fastkernel.h"
