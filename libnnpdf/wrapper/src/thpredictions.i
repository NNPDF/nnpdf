%module(package="NNPDF") thpredictions
 %{
#define SWIG_FILE_WITH_INIT
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/thpredictions.h"
 %}

/*This must be included BEFORE anyrhing from the lib*/
%include "include/numpy.i"
%init %{
  import_array();
%}
%include "std_string.i" 
%include "std_vector.i" 
%include "common.i"


%include "fkgenerator.i"
%include "pdfset.i"

%include "experiments.i"

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"
%include "include/real_typemap.i"

%apply (NNPDF::real** ARGOUTVIEWM_ARRAY1, int* DIM1) {(NNPDF::real** data, int* n)}
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2){(NNPDF::real** data, int*m, int* n)}
%include "../../src/NNPDF/thpredictions.h"

%extend NNPDF::ThPredictions{

void get_data(NNPDF::real **data, int* m, int* n){
    *m = $self->GetNData();
    *n = $self->GetNPdf();
    int len = (*m) * (*n);
    NNPDF::real * result = new NNPDF::real[len];
    NNPDF::real* obs = $self->GetObs();
    std::copy(obs, obs+len ,result);
    *data = result;
}

void get_cv (NNPDF::real  **data, int* n){
    int len = $self->GetNData();
    NNPDF::real * result = new NNPDF::real[len];
    for (int i = 0; i < len; i++){
        result[i] = self->GetObsCV(i);
    }
    *data = result;
    *n = len;
}

void get_error (NNPDF::real  **data, int* n){
    int len = $self->GetNData();
    NNPDF::real * result = new NNPDF::real[len];
    for (int i = 0; i < len; i++){
        result[i] = self->GetObsError(i);
    }
    *data = result;
    *n = len;
}

%pythoncode{

def __len__(self):
    return self.GetNData()


def __iter__(self):
    for i in range(len(self)):
        yield self.GetObsCV(i), self.GetObsError(i)

}
}
