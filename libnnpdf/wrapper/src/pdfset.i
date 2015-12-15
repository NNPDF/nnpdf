%module(package="NNPDF") pdfset
 %{
#include "../../src/NNPDF/pdfset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/pdfset.h"
