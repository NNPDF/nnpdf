%module(package="NNPDF") thpredictions
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/thpredictions.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "fkgenerator.i"
%import "pdfset.i"

%import "experiments.i"

/* Parse the header file to generate wrappers */

%include "include/excepthandler.i"

%include "../../src/NNPDF/thpredictions.h"
