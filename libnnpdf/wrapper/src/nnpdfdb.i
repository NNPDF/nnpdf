%module(package="NNPDF") nnpdfdb
 %{
#include "NNPDF/exceptions.h"
#include "NNPDF/nnpdfdb.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "std_map.i"
%include "common.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"
%template(map_string_string) std::map<std::string,std::string>;

%include "NNPDF/nnpdfdb.h"
