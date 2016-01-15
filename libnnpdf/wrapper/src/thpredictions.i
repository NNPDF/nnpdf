%module(package="NNPDF") thpredictions
 %{
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/thpredictions.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "fkgenerator.i"
%import "pdfset.i"

%import "experiments.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "../../src/NNPDF/thpredictions.h"

%extend NNPDF::ThPredictions{

%pythoncode{

def __len__(self):
    return self.GetNData()


def __iter__(self):
    for i in range(len(self)):
        yield self.GetObsCV(i), self.GetObsError(i)

}
}
