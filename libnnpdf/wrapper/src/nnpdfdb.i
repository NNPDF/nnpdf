%module(package="NNPDF") nnpdfdb
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/nnpdfdb.h"
 %}

%include "std_string.i" 
%include "std_map.i"
%include "std_vector.i"
%include "common.i"
%include "include/excepthandler.i"

/* Parse the header file to generate wrappers */
 //%template(map_string_string) std::map<std::string, std::string>;

%template(vectorString) std::vector<std::string>;
%include "../../src/NNPDF/nnpdfdb.h"
