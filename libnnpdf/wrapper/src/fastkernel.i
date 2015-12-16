%module(package="NNPDF") fastkernel
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/fastkernel.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "exception.i"

%include "common.i"

/* Parse the header file to generate wrappers */
%template(_string_list) std::vector< std::string >;
%include "../../src/NNPDF/exceptions.h"

int NNPDF::ii = 25;

%exception {
  try {
    $action 
  } catch(NNPDF::RuntimeException &_e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(_e.what()));
  }
}

%include "../../src/NNPDF/fastkernel.h"
