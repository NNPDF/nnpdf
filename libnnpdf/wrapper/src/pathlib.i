%module(package="NNPDF") pathlib
%{
#include "NNPDF/exceptions.h"
#include "NNPDF/pathlib.h"
%}

%include "std_string.i"
%include "common.i"

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/pathlib.h"
