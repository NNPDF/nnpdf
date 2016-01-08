%module(package="NNPDF") nnpdfdb
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/nnpdfdb.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "std_map.i"
%include "common.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"
%template(map_string_string) std::map<std::string,std::string>;
%template(vector_string) std::vector<std::string>;

%include "../../src/NNPDF/nnpdfdb.h"
