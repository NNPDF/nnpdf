%module commondata
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

#include "NNPDF/dataset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%import "commondata.i"
%import "fkset.i"
/* Parse the header file to generate wrappers */
%include "../NNPDF/dataset.h"
