%module(package="NNPDF") lhapdfset
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/lhapdfset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "pdfset.i"
%include "common.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/lhapdfset.h"
