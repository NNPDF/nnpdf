%module(package="NNPDF") thpredictions
 %{
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "../../src/NNPDF/thpredictions.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "fkgenerator.i"
%import "pdfset.i"

%import "experiments.i"

/* Parse the header file to generate wrappers */
%include "../../src/NNPDF/thpredictions.h"
