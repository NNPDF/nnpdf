%module(package="NNPDF") experiments
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/experiments.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%import "dataset.i"
%import "pdfset.i"

/* Parse the header file to generate wrappers */

%template (vector_dataset) std::vector<NNPDF::DataSet>;
%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/experiments.h"
