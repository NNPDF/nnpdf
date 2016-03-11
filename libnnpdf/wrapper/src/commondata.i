%module(package="NNPDF") commondata
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/commondata.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "std_map.i" 
%include "common.i"
%template(map_str_vector_str) std::map<std::string, std::vector<std::string> >;
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/commondata.h"
