%module(package="NNPDF") fkset
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/fkset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"

/* Parse the header file to generate wrappers */

%template (vector_fktable_p) std::vector<NNPDF::FKTable*>;
%feature("autodoc", "3");

%include "include/excepthandler.i"
%include "fkgenerator.i"
%include "../../src/NNPDF/fkset.h"
