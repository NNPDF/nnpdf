%module(package="NNPDF") fkset
 %{

#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/fkset.h"
 %}

%include "exception.i"
%include "std_string.i" 
%include "std_vector.i" 

/* Parse the header file to generate wrappers */

%include "../../src/NNPDF/exceptions.h"
%include "fkgenerator.i"
%include "../../src/NNPDF/fkset.h"
