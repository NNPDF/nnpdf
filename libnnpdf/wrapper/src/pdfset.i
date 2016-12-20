%module(package="NNPDF") pdfset
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/pdfset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/pdfset.h"
