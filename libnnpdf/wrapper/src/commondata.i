%module(package="NNPDF") commondata
 %{
#include "../../src/NNPDF/commondata.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"
%ignore NNPDF::dataInfo();
%ignore NNPDF::dataInfoRaw();

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/commondata.h"
